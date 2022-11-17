package mascot.ode;


import org.apache.commons.math3.util.FastMath;

/**
 * @author Nicola Felix Mueller
 */
public class Euler2ndOrderTransitions  {
    
    
	double epsilon;
	double max_step;
	
	double[] migration_rates;
	int n, n2;
	int[] indicators;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
	double[] sumStates;
	boolean hasIndicators;


	
	public Euler2ndOrderTransitions(double[] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        n = (int)(Math.sqrt(migration_rates.length) + 0.5);
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
    	hasIndicators = false;
	}
	
	public Euler2ndOrderTransitions(double[] migration_rates, int[] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        n = (int)(Math.sqrt(migration_rates.length) + 0.5);
        this.indicators = indicators;
        n2 = indicators.length / 2;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];    	
    	hasIndicators = true;
	}

	
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot) throws Exception{
		if (false){
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
	
	private double updateP (double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot) throws Exception{
		double max_dotdotdot = 0.0;
		for (int i = 0; i < p.length; i++){
			if (FastMath.abs(pDotDotDot[i]) > max_dotdotdot)
				max_dotdotdot = FastMath.abs(pDotDotDot[i]);
		}
		
		double timeStep = FastMath.min(FastMath.pow((epsilon*6/max_dotdotdot), 1.0/3), FastMath.min(duration, max_step));

		double timeStepSquare = timeStep*timeStep*0.5;
		for (int i = 0; i < p.length; i++){
			double new_val = p[i] + pDot[i]*timeStep + pDotDot[i]*timeStepSquare;
			double diff = FastMath.abs(new_val - p[i]);
			while (new_val > 1 || new_val < 0 || diff>0.2){
				timeStep *= 0.9;
				if (timeStep<1e-32){
					System.out.println("values are " + p[i] + " " + pDot[i]  + " " + pDotDot[i] + " " + pDotDotDot[i]);
					System.out.println("index is " + i + " of " + p.length + " " + lineages*states);
					System.out.println("minimum step size of 1e-32 reached");
					throw new Exception("error in the calculation of transition probabilities");
				}
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
		for (int i = 0; i < p.length/states; i ++){
			double linSum = 0;
			for (int j = 0; j < states; j++){
				linSum += p[states*i+j];
			}
			for (int j = 0; j < states; j++){
				p[states*i+j] /= linSum;
			}
		}
	}


    
    public void computeDerivatives(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot) {
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
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j < states; j++){
    			double pj = p[currlin+j];
    			for (int k = j+1; k < states; k++){    				
					// the probability of lineage i being in state j is p[i*nr_states +j]
					migrates = p[currlin+k]*migration_rates[k * n + j] -
							pj*migration_rates[j * n + k];
					pDot[currlin+j] += migrates;
					pDot[currlin+k] -= migrates;
    				
    			}// j    			 
    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[currlin+j] = coal;
    			pDotDotDot[currlin+j] = coal;
    			pDot[currlin+j] +=	p[currlin+j] * coal;
    			
        		
        		for (int k = 0; k < states; k++){
        			double sumCoaltransition = 0.0;
        			for (int l = 0; l < states; l++)
        				sumCoaltransition += tCR[l]*p[lineages*states + i*states*states + k*states + l];
        				
        			double coalTransition = sumCoaltransition - tCR[j];
        			pDot[lineages*states + i*states*states + k*states + j] =
      					 p[lineages*states + i*states*states + k*states + j] * coalTransition;
        			pDotDot[lineages*states + i*states*states + k*states + j] = coalTransition;
        			pDotDotDot[lineages*states + i*states*states + k*states + j] = coalTransition;
        		}
    		}// j
    	}// lineages   
    	
    	// Update the transition p dots
    	for (int i = 0; i < lineages; i++){
	    	for (int j = 0; j < states; j++){ // initial value
	    		for (int k = 0; k < states; k++){ // which current state
	    			double pk = p[lineages*states + i*states*states + j*states + k];
	    			for (int l = k+1; l < states; l++){  // which flow to state
	    				 // which lineages
    					migrates = p[lineages*states + i*states*states + j*states + l]*migration_rates[l * n + k] -
    							pk*migration_rates[k * n + l];
    					
    					pDot[lineages*states + i*states*states + j*states + k] += migrates;
						pDot[lineages*states + i*states*states + j*states + l] -= migrates;
	    			}// l
	    		}// k
	    	}// j
    	}
    }   
    
    public void computeSecondDerivate(double[] p, double[] pDot, double[] pDotDot){    	
    	
		double migrates;
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
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			double pj = pDot[currlin+j];
    			for (int k = j+1; k<states; k++){
					// the probability of lineage i being in state j is p[i*nr_states +j]
					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
							pj*migration_rates[j * n + k];
					pDotDot[currlin+j] += migrates;
					pDotDot[currlin+k] -= migrates;

    			}// j    			
    			// Calculate the Derivate of p:
    			pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    			
        		for (int k = 0; k < states; k++){ 	
            		double pCoalRateTransition = 0.0;    		
            		for (int l = 0; l<states; l++)
            			pCoalRateTransition += 
            				coalescent_rates[l] * 
            				(pDot[lineages*states + i*states*states + k*states + l] * (sumStates[l] - p[currlin+l]) 
            						+ p[lineages*states + i*states*states + k*states + l]* (sumDotStates[l] - pDot[currlin+j]));

        			pDotDot[lineages*states + i*states*states + k*states + j] *=
      					pDot[lineages*states + i*states*states + k*states + j];
        			pDotDot[lineages*states + i*states*states + k*states + j] += 
        				p[lineages*states + i*states*states + k*states + j]
        						*(pCoalRateTransition - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
        		}

    		}// j    		
    	}// lineages  
    	
    	
    	// Update the transition p dots
    	for (int i = 0; i < lineages; i++){
	    	for (int j = 0; j < states; j++){ // initial value
	    		for (int k = 0; k < states; k++){ // which current state
	    			double pk = pDot[lineages*states + i*states*states + j*states + k];
	    			for (int l = k+1; l < states; l++){  // which flow to state
	    				 // which lineages
    					migrates = pDot[lineages*states + i*states*states + j*states + l]*migration_rates[l * n + k] -
    							pk*migration_rates[k * n + l];
    					
    					pDotDot[lineages*states + i*states*states + j*states + k] += migrates;
    					pDotDot[lineages*states + i*states*states + j*states + l] -= migrates;
	    			}// l
	    		}// k
	    	}// j
    	}
    }
    
    public void approximateThirdDerivate(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
    	double migrates;
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		// Calculate the probability of a lineage changing states
    		int currlin = states*i;
    		for (int j = 0; j<states; j++){
    			pDotDotDot[currlin+j] *= pDotDot[currlin+j];
        		for (int k = 0; k < states; k++){ 			
        			pDotDotDot[lineages*states + i*states*states + j*states + k] *=
        					pDotDot[lineages*states + i*states*states + j*states + k];
        		}
    		}
    		
    		for (int j = 0; j<states; j++){
    			double pj = pDotDot[currlin+j];
    			for (int k = j+1; k<states; k++){
					migrates = pDotDot[currlin+k]*migration_rates[k * n + j] -
							pj*migration_rates[j * n + k];
					pDotDotDot[currlin+j] += migrates;
					pDotDotDot[currlin+k] -= migrates;
    				
    			}// j      			
    		}// j    		
    	}// lineages  
    	
    	// Update the transition p dots
    	for (int i = 0; i < lineages; i++){
	    	for (int j = 0; j < states; j++){ // initial value
	    		for (int k = 0; k < states; k++){ // which current state
	    			double pk = pDotDot[lineages*states + i*states*states + j*states + k];
	    			for (int l = k+1; l < states; l++){  // which flow to state
	    				 // which lineages
    					migrates = pDotDot[lineages*states + i*states*states + j*states + l]*migration_rates[l * n + k] -
    							pk*migration_rates[k * n + l];
    					
    					pDotDotDot[lineages*states + i*states*states + j*states + k] += migrates;
    					pDotDotDot[lineages*states + i*states*states + j*states + l] -= migrates;
	    			}// l
	    		}// k
	    	}// j
    	}

    }


    
    public void computeDerivativesIndicators(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot) {
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
    		
    		
    		for (int j = 0; j < indicators.length; j++){
    			migrates = p[currlin+indicators[j*n2+0]]*migration_rates[indicators[j*n2+0] * n + indicators[j*n2+1]];
    			pDot[currlin+indicators[j*n2+1]] += migrates;
    			pDot[currlin+indicators[j*n2+0]] -= migrates;  			
    		}
    		
    		for (int j = 0; j < states; j++){    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[currlin+j] = coal;
    			pDotDotDot[currlin+j] = coal;
    			pDot[currlin+j] +=	p[currlin+j] * coal;
    			
        		
        		for (int k = 0; k < states; k++){
        			double sumCoaltransition = 0.0;
        			for (int l = 0; l < states; l++)
        				sumCoaltransition += tCR[l]*p[lineages*states + i*states*states + k*states + l];
        				
        			double coalTransition = sumCoaltransition - tCR[j];
        			pDot[lineages*states + i*states*states + k*states + j] =
      					 p[lineages*states + i*states*states + k*states + j] * coalTransition;
        			pDotDot[lineages*states + i*states*states + k*states + j] = coalTransition;
        			pDotDotDot[lineages*states + i*states*states + k*states + j] = coalTransition;
        		}
    		}// j
    	}// lineages   
    	
    	// Update the transition p dots
    	for (int i = 0; i < lineages; i++){
	    	for (int j = 0; j < states; j++){ // initial value
	    		for (int k = 0; k < indicators.length; k++){
	    			migrates 
	    				= p[lineages*states + i*states*states + j*states+indicators[k*n2+0]]
	    						*migration_rates[indicators[k*n2+0] * n + indicators[j*n2+1]];
	    			pDot[lineages*states + i*states*states + j*states+indicators[k*n2+1]] += migrates;
	    			pDot[lineages*states + i*states*states + j*states+indicators[k*n2+0]] -= migrates;  			
	    		}
	    	}// j
    	}
    }   
    
    public void computeSecondDerivateIndicators(double[] p, double[] pDot, double[] pDotDot){    	
    	
		double migrates;
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
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j < indicators.length; j++){
    			migrates = pDot[currlin+indicators[j*n2+0]]*migration_rates[indicators[j*n2+0] * n + indicators[j*n2+1]];
    			pDotDot[currlin+indicators[j*n2+1]] += migrates;
    			pDotDot[currlin+indicators[j*n2+0]] -= migrates;  			
    		}

    		    		
    		for (int j = 0; j<states; j++){
    			// Calculate the Derivate of p:
    			pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    			
        		for (int k = 0; k < states; k++){ 	
            		double pCoalRateTransition = 0.0;    		
            		for (int l = 0; l<states; l++)
            			pCoalRateTransition += 
            				coalescent_rates[l] * 
            				(pDot[lineages*states + i*states*states + k*states + l] * (sumStates[l] - p[currlin+l]) 
            						+ p[lineages*states + i*states*states + k*states + l]* (sumDotStates[l] - pDot[currlin+j]));

//	    			if ((lineages*states + i*states*states + k*states + j)==5652)
//	    				System.out.println("1 " + pDotDot[lineages*states + i*states*states + k*states + j]);

        			pDotDot[lineages*states + i*states*states + k*states + j] *=
      					pDot[lineages*states + i*states*states + k*states + j];
//	    			if ((lineages*states + i*states*states + k*states + j)==5652)
//	    				System.out.println("2 " + pDotDot[lineages*states + i*states*states + k*states + j]);
	    			pDotDot[lineages*states + i*states*states + k*states + j] += 
        				p[lineages*states + i*states*states + k*states + j]
        						*(pCoalRateTransition - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
//	    			if ((lineages*states + i*states*states + k*states + j)==5652)
//	    				System.out.println("3 " + pDotDot[lineages*states + i*states*states + k*states + j]);
       		}

    		}// j    		
    	}// lineages 
//		if (pDot.length>5652)
//			System.out.println("4 " + pDotDot[5652] + " " +  pDot[5652]);
   	
    	// Update the transition p dots
    	for (int i = 0; i < lineages; i++){
	    	for (int j = 0; j < states; j++){ // initial value
	    		for (int k = 0; k < indicators.length; k++){
	    			double psource = pDot[lineages*states + i*states*states + j*states+indicators[k*n2+0]];
	    			if (psource>0.0){
		    			migrates = psource*migration_rates[indicators[k*n2+0] * n + indicators[k*n2+1]];
//		    			if ((lineages*states + i*states*states + j*states+indicators[k][0])==5652)
//		    				System.out.println("mig + " + migrates + " " + pDotDot[lineages*states + i*states*states + j*states+indicators[k][0]] + " " + pDot[lineages*states + i*states*states + j*states+indicators[k][0]]);
//		    			if ((lineages*states + i*states*states + j*states+indicators[k][1])==5652)
//		    				System.out.println("mig - " + migrates + " " + pDotDot[lineages*states + i*states*states + j*states+indicators[k][1]] + " " + pDot[lineages*states + i*states*states + j*states+indicators[k][0]]);
		    				
		    			pDotDot[lineages*states + i*states*states + j*states+indicators[k*n2+1]] += migrates;
		    			pDotDot[lineages*states + i*states*states + j*states+indicators[k*n2+0]] -= migrates;  
	    			}
	    		}
	    	}// j
    	}

    	
    }
    
    public void approximateThirdDerivateIndicators(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
    	double migrates;
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		// Calculate the probability of a lineage changing states
    		int currlin = states*i;
    		for (int j = 0; j<states; j++){
    			pDotDotDot[currlin+j] *= pDotDot[currlin+j];
        		for (int k = 0; k < states; k++){ 			
        			pDotDotDot[lineages*states + i*states*states + j*states + k] *=
        					pDotDot[lineages*states + i*states*states + j*states + k];
        		}
    		}
    		for (int j = 0; j < indicators.length; j++){
    			migrates = pDotDot[currlin+indicators[j*n2+0]]*migration_rates[indicators[j*n2+0] * n + indicators[j*n2+1]];
    			pDotDotDot[currlin+indicators[j*n2+1]] += migrates;
    			pDotDotDot[currlin+indicators[j*n2+0]] -= migrates;  			
    		}

       	}// lineages  
    	
    	// Update the transition p dots
    	for (int i = 0; i < lineages; i++){
	    	for (int j = 0; j < states; j++){ // initial value
	    		for (int k = 0; k < indicators.length; k++){
	    			migrates = pDotDot[lineages*states + i*states*states + j*states+indicators[k*n2+0]]*migration_rates[indicators[k*n2+0] * n + indicators[k*n2+1]];
	    			pDotDotDot[lineages*states + i*states*states + j*states+indicators[k*n2+1]] += migrates;
	    			pDotDotDot[lineages*states + i*states*states + j*states+indicators[k*n2+0]] -= migrates;  			
	    		}
	    	}// j
    	}
    }

  
}