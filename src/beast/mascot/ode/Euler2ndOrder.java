package beast.mascot.ode;


import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

public class Euler2ndOrder {

	double epsilon;
	double max_step;
	
	double[][] migration_rates;
	int[][] indicators;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
	double[] sumStates;
	boolean hasIndicators;

	
	public Euler2ndOrder(double[][] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.coalescent_rates = new double[states];
        for (int i = 0; i < states; i++){
        	this.coalescent_rates[i] = 2*coalescent_rates[i];
        }
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
    	hasIndicators = false;
	}
	
	public Euler2ndOrder(double[][] migration_rates, int[][] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.indicators = indicators;
    	this.coalescent_rates = coalescent_rates;       
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];
    	hasIndicators = true;
	}

	
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		if (hasIndicators){
			while (duration > 0){
				pDot = new double[pDot.length];
				pDotDot[pDot.length-1] = 0;
				computeDerivativesIndicators(p, pDot, pDotDot, pDotDotDot);
				computeSecondDerivateIndicators(p, pDot, pDotDot);
				approximateThirdDerivateIndicators(p, pDot, pDotDot, pDotDotDot);
				duration = updateP(duration, p,  pDot, pDotDot, pDotDotDot);
			}
		}else{
			while (duration > 0){
				pDot = new double[pDot.length];
				pDotDot[pDot.length-1] = 0;
				computeDerivatives(p, pDot, pDotDot, pDotDotDot);
				computeSecondDerivate(p, pDot, pDotDot);
				approximateThirdDerivate(p, pDot, pDotDot, pDotDotDot);
				duration = updateP(duration, p,  pDot, pDotDot, pDotDotDot);
			}			
		}
	}	
	
	int count = 1;
	private double updateP (double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		double max_dotdotdot = 0.0;
		for (int i = 0; i < (p.length-1); i++){
			if (FastMath.abs(pDotDotDot[i]) > max_dotdotdot)
				max_dotdotdot = FastMath.abs(pDotDotDot[i]);
		}
				
		double timeStep = FastMath.min(FastMath.pow((epsilon*6/max_dotdotdot), 1.0/3), FastMath.min(duration, max_step));

		double timeStepSquare = timeStep*timeStep*0.5;
		for (int i = 0; i < (p.length-1); i++){
			double new_val = p[i] + pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
			double diff = FastMath.abs(new_val - p[i]);
			while (new_val > 1 || new_val < 0 || diff>0.2){
				timeStep *= 0.9;
				timeStepSquare = timeStep*timeStep*0.5;
				new_val = p[i] + pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
				diff = FastMath.abs(new_val - p[i]);
			}			
		}
		doUpdating(timeStep, timeStepSquare, p, pDot, pDotDot);
		duration -= timeStep;
		return duration;
		
		
	}
	
	private void doUpdating(double timeStep, double timeStepSquare, double[] p, double[] pDot, double[] pDotDot){
		for (int i = 0; i < p.length; i++)
			p[i] += pDot[i]*timeStep + pDotDot[i]*timeStepSquare;	
		
		// normalize to ensure stability
		for (int i = 0; i < lineages; i ++){
			double linSum = 0;
			for (int j = 0; j < states; j++){
				linSum += p[states*i+j];
			}
			for (int j = 0; j < states; j++){
				p[states*i+j] /= linSum;
			}
		}
	}
	
    
	public void computeDerivatives (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot) {
    	double migrates;
    	// Compute the sum of line state probabilities for each state
     	sumStates = new double[states];
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
				sumStates[j] += p[states*i+j];    			
    		
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){
    		double[] tCR =  new double[states];
    		double sumCoal = 0;
    		int currlin = states*i;
    		for (int j = 0; j<states; j++){
    			tCR[j] = coalescent_rates[j] *  (sumStates[j] - p[currlin+j]);
    			sumCoal += p[currlin+j]*tCR[j];
    		}
    		pDot[pDot.length-1] -= sumCoal;
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j < states; j++){
    			double pj = p[currlin+j];
    			for (int k = j+1; k < states; k++){    
    				
					// the probability of lineage i being in state j is p[i*nr_states +j]
					migrates = p[currlin+k]*migration_rates[k][j] -
							pj*migration_rates[j][k];
					pDot[currlin+j] += migrates;
					pDot[currlin+k] -= migrates;
//					if (j==0 && k==4)
//						System.out.println(migrates);
    				
    			}// j    			 
    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[currlin+j] = coal;
    			pDotDotDot[currlin+j] = coal;
    			pDot[currlin+j] +=	p[currlin+j] * coal;
    		}// j
    	}// lineages    
    	
		pDot[pDot.length-1]  /= 2;
    }
        
    public void computeSecondDerivate(double[] p, double[] pDot, double[] pDotDot){    	
    	
    	double[] sumDotStates = new double[states];
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
    			sumDotStates[j] += pDot[states*i+j];      	
	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));     
    		
    		for (int j = 0; j<states; j++)
    			pDotDot[currlin+j] *= pDot[currlin+j];
    		
    		pDotDot[pDot.length-1] -= pCoalRate;
    		
    		double migrates;
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			double pj = pDot[currlin+j];
    			for (int k = j+1; k<states; k++){
					// the probability of lineage i being in state j is p[i*nr_states +j]
					migrates = pDot[currlin+k]*migration_rates[k][j] -
							pj*migration_rates[j][k];
					pDotDot[currlin+j] += migrates;
					pDotDot[currlin+k] -= migrates;
    			}// j    			
    			// Calculate the Derivate of p:
    			pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    		}// j    		
    	}// lineages    	
		pDotDot[pDot.length-1] /= 2;		

    }
    
    public void approximateThirdDerivate(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
    	double migrates;
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		// Calculate the probability of a lineage changing states
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pDotDotDot[currlin+j] *= pDotDot[currlin+j];
    		
    		for (int j = 0; j<states; j++){
    			double pj = pDotDot[currlin+j];
    			for (int k = j+1; k<states; k++){
					migrates = pDotDot[currlin+k]*migration_rates[k][j] -
							pj*migration_rates[j][k];
					pDotDotDot[currlin+j] += migrates;
					pDotDotDot[currlin+k] -= migrates;
    				
    			}// j      			
    		}// j    		
    	}// lineages        	
    }
   
	public void computeDerivativesIndicators (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot) {
    	double migrates;
    	// Compute the sum of line state probabilities for each state
     	sumStates = new double[states];
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
				sumStates[j] += p[states*i+j];    			
    		
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){
    		double[] tCR =  new double[states];
    		double sumCoal = 0;
    		int currlin = states*i;
    		for (int j = 0; j<states; j++){
    			tCR[j] = coalescent_rates[j] *  (sumStates[j] - p[currlin+j]);
    			sumCoal += p[currlin+j]*tCR[j];
    		}
    		pDot[pDot.length-1] -= sumCoal;
    		for (int j = 0; j < states; j++){    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[currlin+j] = coal;
    			pDotDotDot[currlin+j] = coal;
    			pDot[currlin+j] +=	p[currlin+j] * coal;
    		}// j

    	}
    		// Calculate the probability of a lineage changing states
		for (int j = 0; j < indicators.length; j++){
			int source = indicators[j][0];
			int sink = indicators[j][1];
			double mrate = migration_rates[source][sink];
	    	for (int i = 0; i<lineages; i++){
				migrates = p[states*i+source]*mrate;
				pDot[states*i+sink] += migrates;
				pDot[states*i+source] -= migrates;  			
	    	}
    	}    
    	
		pDot[pDot.length-1]  /= 2;
    }
        
    public void computeSecondDerivateIndicators(double[] p, double[] pDot, double[] pDotDot){    	
    	
    	double[] sumDotStates = new double[states];
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
    			sumDotStates[j] += pDot[states*i+j];      	
	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));     
    		
    		for (int j = 0; j<states; j++)
    			pDotDot[currlin+j] *= pDot[currlin+j];
    		
    		pDotDot[pDot.length-1] -= pCoalRate;
    		

    		for (int j = 0; j<states; j++){
    			pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    		}// j    		
    	}// lineages 
    	
		double migrates;
		
		
		// Calculate the probability of a lineage changing states
		for (int j = 0; j < indicators.length; j++){
			int source = indicators[j][0];
			int sink = indicators[j][1];
			double mrate = migration_rates[source][sink];
	    	for (int i = 0; i<lineages; i++){
				migrates = pDot[states*i+source]*mrate;
				pDotDot[states*i+sink] += migrates;
				pDotDot[states*i+source] -= migrates;  			
	    	}
    	}    
		pDotDot[pDot.length-1] /= 2;		

    }
    
    public void approximateThirdDerivateIndicators(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
    	double migrates;
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		// Calculate the probability of a lineage changing states
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pDotDotDot[currlin+j] *= pDotDot[currlin+j];
    	}
    	
		// Calculate the probability of a lineage changing states
		for (int j = 0; j < indicators.length; j++){
			int source = indicators[j][0];
			int sink = indicators[j][1];
			double mrate = migration_rates[source][sink];
	    	for (int i = 0; i<lineages; i++){
				migrates = pDotDot[states*i+source]*mrate;
				pDotDotDot[states*i+sink] += migrates;
				pDotDotDot[states*i+source] -= migrates;  			
	    	}
    	}    
    }

}


