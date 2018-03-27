package beast.mascot.ode;


import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

public class Euler2ndOrderMem {

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

	
	public Euler2ndOrderMem(double[][] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
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
	
	public Euler2ndOrderMem(double[][] migration_rates, int[][] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
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

	public Euler2ndOrderMem(int[] multiplicator, double[][] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
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
	
	public Euler2ndOrderMem(int[] multiplicator, double[][] migration_rates, int[][] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
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
			computeDerivatives(p, pDot, pDotDot, pDotDotDot);
			computeSecondDerivate(p, pDot, pDotDot);
			approximateThirdDerivate(p, pDot, pDotDot, pDotDotDot);
			duration = updateP(duration, p,  pDot, pDotDot, pDotDotDot);			
		}			
	}	
	
	private void getDerivatives(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		
		get1(p, pDot, pDotDot, pDotDotDot);
		get2(p, pDot, pDotDot, pDotDotDot);
		get3(p, pDot, pDotDot, pDotDotDot);
				
	}
	private void get1(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		computeDerivatives(p, pDot, pDotDot, pDotDotDot);
	}
	private void get2(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		computeSecondDerivate(p, pDot, pDotDot);
	}
	private void get3(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		approximateThirdDerivate(p, pDot, pDotDot, pDotDotDot);
	}
	
	private double updateP (double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		double max_dotdotdot = 0.0;
		for (int i = 0; i < dimension; i++){
			max_dotdotdot = FastMath.max(max_dotdotdot, FastMath.abs(pDotDotDot[i]));
		}
		
				
		double timeStep = FastMath.min(FastMath.pow((epsilon*6/max_dotdotdot), 1.0/3), FastMath.min(duration, max_step));

		iterations=0;
		double timeStepSquare = timeStep*timeStep*0.5;
		
		double new_val;

		for (int i = 0; i < dimension; i++){
			new_val = p[i] + pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
			while (new_val > 1 || new_val < 0){
				iterations++;
				timeStep *= 0.5;
				timeStepSquare = timeStep*timeStep*0.5;
				new_val = p[i] + pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
				if (iterations>100){
					System.err.println("too many iterations, return negative infinity");
					p[dimension] = Double.NEGATIVE_INFINITY;
					return 0.0;
				}
			}			
		}
		for (int i = 0; i < dimension; i++){
			p[i] += pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
		}	
		p[dimension] = p[dimension] + pDot[dimension]*timeStep + pDotDot[dimension]*timeStepSquare;
		
		duration -= timeStep;
		return duration;		
	}
	
	private void doUpdating(double timeStep, double timeStepSquare, double[] p, double[] pDot, double[] pDotDot){
		for (int i = 0; i < dimension+1; i++)
			p[i] += pDot[i]*timeStep + pDotDot[i]*timeStepSquare;	
		
//		// normalize to ensure stability
//		double linSum;
//		for (int i = 0; i < lineages; i ++){
//			linSum = 0;
//			for (int j = 0; j < states; j++){
//				if (p[states*i+j]>=0.0){
//					linSum += p[states*i+j];
//				}else{
//					if (p[states*i+j]>=-1e-300){
//						p[states*i+j] = 0.0;
//					}else{
//						System.err.println("value below zero " + p[states*i+j]);
//						p[dimension] = Double.NEGATIVE_INFINITY;
//						return;
//					}
//				}
//			}
//			for (int j = 0; j < states; j++){
//				p[states*i+j] /= linSum;
//			}
//		}
	}
	    
	private void computeDerivatives (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot) {
    	// Compute the sum of line state probabilities for each state
     	sumStates = new double[states];
     	computeSumStates(p, sumStates);
    		
    	// Caluclate the change in the lineage state probabilities for every lineage in every state     	
     	coalFirst(p, pDot, pDotDot, pDotDotDot);    	
     	migrate(p, pDot);    	    	
		pDot[dimension]  /= 2;
    }
	
	private void computeSecondDerivate (double[] p, double[] pDot, double[] pDotDot){  
    	double[] sumDotStates = new double[states];
    	computeSumStates(pDot, sumDotStates);
    	coalSecond(p, pDot, pDotDot, sumDotStates);     	
    	migrate(pDot, pDotDot);
		pDotDot[dimension] /= 2;
    }
	
    private void approximateThirdDerivate (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		// Calculate the probability of a lineage changing states
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pDotDotDot[currlin+j] *= pDotDot[currlin+j];
    	}   
    	
    	migrate(pDotDot, pDotDotDot);
    }
    
	private void coalFirst(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
//		for (int i = 0; i < dimension; i++){
//			pDot[i] = 0;
//			pDotDot[i] = 0;
//			pDotDotDot[i] = 0;
//		}
     	double[] tCR =  new double[states];
     	double sumCoal,coal;
     	int currlin;
     	pDot[dimension] = 0;
     	pDotDot[dimension] = 0;

     	if (hasMultiplicator){
	    	for (int i = 0; i<lineages; i++){    		
	    		sumCoal = 0;
	    		currlin = states*i;
	    		for (int j = 0; j<states; j++){
	    			tCR[j] = coalescent_rates[j] *  (sumStates[j] - p[currlin+j]);
	    			sumCoal += p[currlin+j]*tCR[j];
	    		}
	    		
	     		pDot[dimension] -= multiplicator[i]*sumCoal;
	         	
	    		for (int j = 0; j < states; j++){    			
	    			// Calculate the Derivate of p:
	    			coal = sumCoal - tCR[j];
	    			pDot[currlin+j] = p[currlin+j] * coal;
	    			pDotDot[currlin+j] = coal;
	    			pDotDotDot[currlin+j] = coal;
	    		}
	    	}  
     	}else{
	    	for (int i = 0; i<lineages; i++){    		
	    		sumCoal = 0;
	    		currlin = states*i;
	    		for (int j = 0; j<states; j++){
	    			tCR[j] = coalescent_rates[j] *  (sumStates[j] - p[currlin+j]);
	    			sumCoal += p[currlin+j]*tCR[j];
	    		}
	    		
	     		pDot[dimension] -= sumCoal;
	         	
	    		for (int j = 0; j < states; j++){    			
	    			// Calculate the Derivate of p:
	    			coal = sumCoal - tCR[j];
	    			pDot[currlin+j] = p[currlin+j] * coal;
	    			pDotDot[currlin+j] = coal;
	    			pDotDotDot[currlin+j] = coal;
	    		}
	    	}  
	    		
	    }

	}
		
	private void coalSecond(double[] p, double[] pDot, double[] pDotDot, double[] sumDotStates){
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	double pCoalRate;
    	int currlin;
    	for (int i = 0; i<lineages; i++){    		
    		pCoalRate = 0.0;    		
    		currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));     
    		
    		for (int j = 0; j<states; j++)
    			pDotDot[currlin+j] *= pDot[currlin+j];
    		
    		if (hasMultiplicator){
    			pDotDot[dimension] -= multiplicator[i]*pCoalRate;
    		}else{
    			pDotDot[dimension] -= pCoalRate;    			
    		}    		

    		for (int j = 0; j<states; j++){
    			pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    		}// j    		
    	}// lineages 

	}
    
    private void migrate(double[] p, double[] pdot){
    	double migrates;
    	if (hasIndicators){
        	int source, sink;
        	double mrate;
			for (int j = 0; j < indicators.length; j++){
				source = indicators[j][0];
				sink = indicators[j][1];
				mrate = migration_rates[source][sink];
		    	for (int i = 0; i<lineages; i++){
					migrates = p[states*i+source]*mrate;
					pdot[states*i+sink] += migrates;
					pdot[states*i+source] -= migrates;  			
		    	}
	    	}    
    	}else{
    		int currlin;
	    	for (int i = 0; i<lineages; i++){
	    		currlin = states*i;
				for (int j = 0; j < states; j++){
					for (int k = 0; k < states; k++){  
						migrates = p[currlin+j]*migration_rates[j][k];
						pdot[currlin+k] += migrates;
						pdot[currlin+j] -= migrates;  					    			 			
					}
	    		}
	    	}
    	}    	
    }

	private void computeSumStates(double[] p, double[] sumStates){
		int currlin;
    	if (hasMultiplicator){    		
	    	for (int i = 0; i<lineages; i++){
	    		currlin = states*i;	    		
	    		for (int j = 0; j<states; j++){
	    			sumStates[j] += multiplicator[i]*p[currlin+j];
	    		}
	    	}
    	}else{
	    	for (int i = 0; i<lineages; i++){
	    		currlin = states*i;
	    		for (int j = 0; j<states; j++){
	    			sumStates[j] += p[currlin+j];
	    		}
	    	}
    	}
		
	}

}


