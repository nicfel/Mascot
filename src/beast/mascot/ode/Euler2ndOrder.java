package beast.mascot.ode;

import java.util.Arrays;

import org.apache.commons.math4.util.FastMath;

public class Euler2ndOrder {

	double epsilon = 0.001;
	
	double[][] migration_rates;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
	double[] sumStates;
	double[] totCoalRate;	

	
	public Euler2ndOrder(double[][] migration_rates, double[] coalescent_rates, int lineages, int states) {
        this.migration_rates = migration_rates;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
	}

	
	
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		while (duration > 0){
			pDot = new double[pDot.length];
			pDotDot = new double[pDot.length];
			pDotDotDot = new double[pDot.length];
			computeDerivatives(p, pDot);
			computeSecondDerivate(p, pDot, pDotDot);
			approximateThirdDerivate(p, pDot, pDotDot, pDotDotDot);
			duration = updateP(duration, p,  pDot, pDotDot, pDotDotDot);
		}
	}	
	
	private double updateP (double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		double max_dotdotdot = 0.0;
		for (int i = 0; i < (p.length-1); i++){
			if (FastMath.abs(pDotDotDot[i]) > max_dotdotdot)
				max_dotdotdot = FastMath.abs(pDotDotDot[i]);
		}
		
		double timeStep = FastMath.min(FastMath.pow((epsilon*6/max_dotdotdot), 1.0/3), duration);

		double timeStepSquare = timeStep*timeStep*0.5;
		for (int i = 0; i < (p.length-1); i++){
			double new_val = p[i] + pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
			
			double lower = p[i]*0.1;
			double upper = (1-p[i])*0.9 + p[i];
			
			while (new_val > upper || new_val < lower){
				timeStep *= 0.9;
				timeStepSquare = timeStep*timeStep*0.5;
				new_val = p[i] + pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
			}			
		}
		doUpdating(timeStep, timeStepSquare, p, pDot, pDotDot);
		duration -= timeStep;
		return duration;
		
		
	}
	
	private void doUpdating(double timeStep, double timeStepSquare, double[] p, double[] pDot, double[] pDotDot){
		for (int i = 0; i < (p.length); i++)
			p[i] += pDot[i]*timeStep + pDotDot[i]*timeStepSquare;		
	}
	
    public void computeDerivatives (double[] p, double[] pDot) {
    	totCoalRate = new double[(p.length-1)/states];
    	double migrates;
    	// Compute the sum of line state probabilities for each state
     	sumStates = new double[states];
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
				sumStates[j] += p[states*i+j];    			
    		
    	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){
    		double tCR = 0.0;
    		for (int j = 0; j<states; j++)
    			tCR += 2*coalescent_rates[j] * p[states*i+j]* (sumStates[j] - p[states*i+j]);     
    		
    		totCoalRate[i] = tCR;
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j < states; j++){
    			for (int k = j+1; k < states; k++){    				
					// the probability of lineage i being in state j is p[i*nr_states +j]
					migrates = p[states*i+k]*migration_rates[k][j] -
							p[states*i+j]*migration_rates[j][k];
					pDot[states*i+j] += migrates;
					pDot[states*i+k] -= migrates;
    				
    			}// j    			 
    			
    			// Calculate the Derivate of p:
    			pDot[states*i+j] +=	p[states*i+j] * (tCR - 2*coalescent_rates[j] * (sumStates[j] - p[states*i+j]));
    		}// j
    	}// lineages    
    	
		pDot[pDot.length-1]  = getTotalCoalescent(sumStates, p);
   	
    }
    
    
    public void computeSecondDerivate(double[] p, double[] pDot, double[] pDotDot){    	
    	
    	double[] sumDotStates = new double[states];
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
    			sumDotStates[j] += pDot[states*i+j];      	
	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		for (int j = 0; j<states; j++)
    			pCoalRate += 2*coalescent_rates[j] * pDot[states*i+j]* (sumStates[j] - p[states*i+j]);     
    		for (int j = 0; j<states; j++)
    			pCoalRate += 2*coalescent_rates[j] * p[states*i+j]* (sumDotStates[j] - pDot[states*i+j]);     
   		
    		
    		double migrates;
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			for (int k = j+1; k<states; k++){
					// the probability of lineage i being in state j is p[i*nr_states +j]
					migrates = pDot[states*i+k]*migration_rates[k][j] -
							pDot[states*i+j]*migration_rates[j][k];
					pDotDot[states*i+j] += migrates;
					pDotDot[states*i+k] -= migrates;

    			}// j    			
    			// Calculate the Derivate of p:
    			pDotDot[states*i+j] += pDot[states*i+j]*totCoalRate[i] + p[states*i+j]*pCoalRate -
    					p[states*i+j] * 2*coalescent_rates[j] * (sumDotStates[j] - pDot[states*i+j])-
						pDot[states*i+j] * 2*coalescent_rates[j] * (sumStates[j] - p[states*i+j]);
    		}// j    		
    	}// lineages    	
		pDotDot[pDot.length-1]  = getSecondCoalescent(sumStates, sumDotStates, p, pDot);		

    }
    
    public void approximateThirdDerivate(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
    	double migrates;
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			for (int k = j+1; k<states; k++){
					migrates = pDotDot[states*i+k]*migration_rates[k][j] -
							pDotDot[states*i+j]*migration_rates[j][k];
					pDotDotDot[states*i+j] += migrates;
					pDotDotDot[states*i+k] -= migrates;
    				
    			}// j      			
    			// Calculate the Derivate of p:
    			pDotDotDot[states*i+j] += pDotDot[states*i+j]* (totCoalRate[i] - 2*coalescent_rates[j] * (sumStates[j] - p[states*i+j]));
    		}// j    		
    	}// lineages        	
    }
   

    private double getTotalCoalescent(double[] sumStates, double[] p_norm){
    	double[] dTdtStates = new double[states];
    	for (int i = 0; i < lineages; i++)
    		for (int j = 0; j<states; j++)
    			dTdtStates[j] += p_norm[i*states+j]*(sumStates[j] - p_norm[i*states+j]);
	
    	
    	double logPdot = 0;
		for (int j = 0; j<states; j++)
			logPdot -= coalescent_rates[j]*dTdtStates[j];
		
		return logPdot;
    }
    
    private double getSecondCoalescent(double[] sumStates, double[] sumDotStates, double[] p_norm, double[] pDot){
    	
    	double[] dTdtStates = new double[states];
    	for (int i = 0; i < lineages; i++)
    		for (int j = 0; j<states; j++)
    			dTdtStates[j] += pDot[i*states+j]*(sumStates[j] - p_norm[i*states+j]);
    	for (int i = 0; i < lineages; i++)
    		for (int j = 0; j<states; j++)
    			dTdtStates[j] += p_norm[i*states+j]*(sumDotStates[j] - pDot[i*states+j]);	
    	
    	double logPdotDot = 0;
		for (int j = 0; j<states; j++)
			logPdotDot -= coalescent_rates[j]*dTdtStates[j];
		
		return logPdotDot;
    }
  


	
}


