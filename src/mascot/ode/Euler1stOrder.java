package mascot.ode;


import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

public class Euler1stOrder {

	double epsilon;
	double max_step;
	
	double[][] migration_rates;
	int[] multiplicator;
	int[][] indicators;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
	double[] sumStates;
	boolean hasIndicators;
	boolean hasMultiplicator;
	
	int iterations;

	
	public Euler1stOrder(double[][] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
    	hasIndicators = false;
    	hasMultiplicator = false;
	}
	
	public Euler1stOrder(double[][] migration_rates, int[][] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
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
    	hasMultiplicator = false;
	}

	public Euler1stOrder(int[] multiplicator, double[][] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.multiplicator = multiplicator;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
    	hasIndicators = false;
    	hasMultiplicator = true; 
	}
	
	public Euler1stOrder(int[] multiplicator, double[][] migration_rates, int[][] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.multiplicator = multiplicator;
        this.indicators = indicators;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
    	hasIndicators = true;
    	hasMultiplicator = true;      	
	}

	
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		while (duration > 0){
	    	
			pDot = new double[pDot.length];
			computeDerivatives(p, pDot, pDotDot, pDotDotDot);
			computeSecondDerivate(p, pDot, pDotDot);
//			approximateThirdDerivate(p, pDot, pDotDot, pDotDotDot);
			duration = updateP(duration, p,  pDot, pDotDot);
			
				
		}			
	}	
	
	private double updateP (double duration, double[] p, double[] pDot, double[] pDotDot){
		double max_dotdotdot = 0.0;
		for (int i = 0; i < (p.length-1); i++){
			if (FastMath.abs(pDotDot[i]) > max_dotdotdot)
				max_dotdotdot = FastMath.abs(pDotDot[i]);
		}
				
		double timeStep = FastMath.min(FastMath.pow((epsilon*6/max_dotdotdot), 1.0/3), FastMath.min(duration, max_step));

		iterations=0;
		for (int i = 0; i < (p.length-1); i++){
			double new_val = p[i] + pDot[i]*timeStep;
			double diff = FastMath.abs(new_val - p[i]);
			while (new_val > 1 || new_val < 0 || diff>0.2){
				iterations++;
				timeStep *= 0.9;
				new_val = p[i] + pDot[i]*timeStep;
				diff = FastMath.abs(new_val - p[i]);
				if (iterations>100){
					System.err.println("too many iterations, return negative infinity");
					p[p.length-1] = Double.NEGATIVE_INFINITY;
					return 0.0;
				}
			}			
		}
		doUpdating(timeStep, p, pDot);
		duration -= timeStep;
		return duration;
		
		
	}
	
	private void doUpdating(double timeStep, double[] p, double[] pDot){
		for (int i = 0; i < p.length; i++)
			p[i] += pDot[i]*timeStep;	
		// normalize to ensure stability
		for (int i = 0; i < lineages; i ++){
			double linSum = 0;
			for (int j = 0; j < states; j++){
				if (p[states*i+j]>=0.0){
					linSum += p[states*i+j];
				}else{
					if (p[states*i+j]>=-1e-300){
						p[states*i+j] = 0.0;
					}else{
						System.err.println("value below zero " + p[states*i+j]);
						p[p.length-1] = Double.NEGATIVE_INFINITY;
						return;
					}
				}
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
     	if (hasMultiplicator){
	    	for (int i = 0; i<lineages; i++)
	    		for (int j = 0; j<states; j++)
					sumStates[j] += multiplicator[i]*p[states*i+j]; 
     	}else{
	    	for (int i = 0; i<lineages; i++)
	    		for (int j = 0; j<states; j++)
					sumStates[j] += p[states*i+j];      		
     	}
    		
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){
    		double[] tCR =  new double[states];
    		double sumCoal = 0;
    		int currlin = states*i;
    		for (int j = 0; j<states; j++){
    			tCR[j] = coalescent_rates[j] *  (sumStates[j] - p[currlin+j]);
    			sumCoal += p[currlin+j]*tCR[j];
    		}
         	if (hasMultiplicator){
         		pDot[pDot.length-1] -= multiplicator[i]*sumCoal;
         	}else{
         		pDot[pDot.length-1] -= sumCoal;        		
         	}
    		for (int j = 0; j < states; j++){    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[currlin+j] = coal;
    			pDotDotDot[currlin+j] = coal;
    			pDot[currlin+j] +=	p[currlin+j] * coal;
    		}// j

    	}
    	
    	
    	// Calculate the probability of a lineage changing states
    	if (hasIndicators){
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
    	}else{
        	for (int i = 0; i<lineages; i++){
        		int currlin = states*i;
        		// Calculate the probability of a lineage changing states
        		for (int j = 0; j < states; j++){
        			double pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k][j] -
    							pj*migration_rates[j][k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j    			 
        		}// j
        	}// lineages       		
    	}
    	
		pDot[pDot.length-1]  /= 2;

    }
        
    public void computeSecondDerivate (double[] p, double[] pDot, double[] pDotDot){  
    	double[] sumDotStates = new double[states];
    	if (hasMultiplicator){
	    	for (int i = 0; i<lineages; i++)
	    		for (int j = 0; j<states; j++)
	    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];   
    	}else{
	    	for (int i = 0; i<lineages; i++)
	    		for (int j = 0; j<states; j++)
	    			sumDotStates[j] += pDot[states*i+j];   
    	}
	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));     
    		
    		for (int j = 0; j<states; j++)
    			pDotDot[currlin+j] *= pDot[currlin+j];
    		
    		if (hasMultiplicator){
    			pDotDot[pDot.length-1] -= multiplicator[i]*pCoalRate;
    		}else{
    			pDotDot[pDot.length-1] -= pCoalRate;    			
    		}    		

    		for (int j = 0; j<states; j++){
    			pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    		}// j    		
    	}// lineages 
    	
		double migrates;
		
		
		// Calculate the probability of a lineage changing states
		if (hasIndicators){
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
		}else{
        	for (int i = 0; i<lineages; i++){
        		int currlin = states*i;
        		// Calculate the probability of a lineage changing states
        		for (int j = 0; j < states; j++){
        			double pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k][j] -
    							pj*migration_rates[j][k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j    			 
        		}// j
        	}// lineages    
			
		}
		pDotDot[pDot.length-1] /= 2;
    }
    
    public void approximateThirdDerivate (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
    	
    	
    	double migrates;
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		// Calculate the probability of a lineage changing states
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pDotDotDot[currlin+j] *= pDotDot[currlin+j];
    	}
    	
		// Calculate the probability of a lineage changing states
    	if (hasIndicators){
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
    	}else{
			for (int j = 0; j < states; j++){
				for (int k = 0; k < states; k++){  
					double mrate = migration_rates[j][k];
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[states*i+j]*mrate;
						pDotDotDot[states*i+k] += migrates;
						pDotDotDot[states*i+j] -= migrates;  			
	    			} 			
				}
	    	}
    	}
    }

}


