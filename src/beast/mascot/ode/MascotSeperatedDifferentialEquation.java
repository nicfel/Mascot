package beast.mascot.ode;


import java.util.Arrays;

import beast.core.Description;


/**
 * @author Nicola Felix Mueller
 */
@Description("Describes the derivates of lineage state probabilities for LISCO as described in Mueller et al., 2016")
public class MascotSeperatedDifferentialEquation  {

	double[][] migration_rates;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
	double[] sumStates;
    
    boolean belowzero = false;

    // constructor
    public MascotSeperatedDifferentialEquation(double[][] migration_rates, double[] coalescent_rates, int lineages , int states){
        this.migration_rates = migration_rates;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        belowzero = false;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];
    }

    public int getDimension() {
        return dimension+1;
    }

    public void computeDerivatives (double[] p, double[] pDot) {
    	// normalize each lineage such that the probability of the lineage not having coalesced
    	// is 1. This is done such that the coalescent rate is calculated conditional on all 
    	// lineages still existing/ not having coalesced
//    	normalizeLineages(p);
    	
    	// Compute the sum of line state probabilities for each state
     	sumStates = new double[states];
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
				sumStates[j] += p[states*i+j];    			
    		
    	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){
    		double totCoalRate = 0.0;
    		for (int j = 0; j<states; j++)
    			totCoalRate += 2*coalescent_rates[j] * p[states*i+j]* (sumStates[j] - p[states*i+j]);     
    		
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			double migrates = 0.0;
    			for (int k = 0; k<states; k++){
    				if (j != k){    					    					
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates += p[states*i+k]*migration_rates[k][j] -
    							p[states*i+j]*migration_rates[j][k];
    				}
    			}// j    			 
    			
    			// Calculate the Derivate of p:
    			pDot[states*i+j] = migrates +
    					p[states*i+j] * (totCoalRate - 2*coalescent_rates[j] * (sumStates[j] - p[states*i+j]));
    		}// j
    	}// lineages    
    	
		pDot[pDot.length-1]  = getTotalCoalescent(sumStates, p);
   	
    }
    
    private void normalizeLineages(double[] p_in){
    	for (int i = 0; i<lineages; i++){
    		double linSum=0.0;
    		for (int j = 0; j<states; j++){
    			linSum += p_in[states*i+j];    			
    		}
    		if (Math.abs(linSum - 1) > 0.000001){
    			System.out.println(Arrays.toString(p_in));
    			System.exit(0);
    		}
    		for (int j = 0; j<states; j++){    			
    			p_in[states*i+j] = p_in[states*i+j]/linSum;   
    		}
    	}
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
    
    public double getLogPDotDot(double[] sumStates, double[] sumDotStates, double[] p_norm, double[] pDot){
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




    public void computeSecondDerivate(double[] p, double[] pDot, double[] pDotDot){
    	double[] sumDotStates = new double[states];
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
    			sumDotStates[j] += pDot[states*i+j];    
    	
		double[] pDoti = new double[states];
		
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){
    		
    		double totPdotCoalRate = 0.0;
    		double totCoalRate = 0.0;
    		for (int j = 0; j<states; j++)
    			totPdotCoalRate += 2*coalescent_rates[j] * p[states*i+j]* (sumStates[j] - p[states*i+j]);     
    		for (int j = 0; j<states; j++)
    			totCoalRate += 2*coalescent_rates[j] * pDot[states*i+j]* (sumStates[j] - p[states*i+j]);     
    		for (int j = 0; j<states; j++)
    			totCoalRate += 2*coalescent_rates[j] * p[states*i+j]* (sumDotStates[j] - pDot[states*i+j]);     
   		
    		
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			double migrates = 0.0;
    			for (int k = 0; k<states; k++){
    				if (j != k){    					    					
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates += pDot[states*i+k]*migration_rates[k][j] -
    							pDot[states*i+j]*migration_rates[j][k];
    				}
    			}// j    			 
    			
    			// Calculate the Derivate of p:
    			pDotDot[states*i+j] = migrates + pDot[states*i+j]*totPdotCoalRate + p[states*i+j]*totCoalRate -
    					p[states*i+j] * 2*coalescent_rates[j] * (sumDotStates[j] - pDot[states*i+j])-
						pDot[states*i+j] * 2*coalescent_rates[j] * (sumStates[j] - p[states*i+j]);
    		}// j
    		
    		
    	}// lineages    
    	
//		pDot[pDot.length-1]  = getTotalCoalescent(sumStates, p);
		pDotDot[pDot.length-1]  = getSecondCoalescent(sumStates, sumDotStates, p, pDot);
		
//		while ( pDotDot[p.length-1]>0){
//			System.out.println( pDotDot[p.length-1]);
//			System.out.println(Arrays.toString(sumStates));
//			System.out.println(Arrays.toString(p));
//			System.out.println(Arrays.toString(sumDotStates));
//			System.out.println(Arrays.toString(pDot));
//			System.exit(0);
//		}

    }
    
    
}