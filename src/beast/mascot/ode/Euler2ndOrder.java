package beast.mascot.ode;


import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

public class Euler2ndOrder implements Euler2ndOrderBase {

	double epsilon;
	double max_step;
	
	double[] migration_rates; // flattened square matrix of migration rates
	int n, n2; // dimension of migration rate matrix and indicators matrix
	int[] multiplicator;
	int[] indicators;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
	double[] sumStates;
	boolean hasIndicators;
	boolean hasMultiplicator;
	double[] tCR;
	double[] sumDotStates;

	int iterations;

	public Euler2ndOrder() {};
	@Override
	public void init(double[] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        n = (int)(Math.sqrt(migration_rates.length) + 0.5);        
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];tCR = new double[states]; sumDotStates = new double[states];   	
    	hasIndicators = false;
    	hasMultiplicator = false;
    	
    	iterations=0;	
	}
	@Override
	public void initWithIndicators(double[] migration_rates, int[] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
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
    	sumStates = new double[states];tCR = new double[states];
    	hasIndicators = true;
    	hasMultiplicator = false;
    	
    	iterations=0;
	}
	
    double [] linProbs_tmpdt;
    double [] linProbs_tmpddt;
    double [] linProbs_tmpdddt;

	@Override
	public void setup(int maxSize) {
		linProbs_tmpdt = new double[maxSize];
		linProbs_tmpddt = new double[maxSize];
		linProbs_tmpdddt = new double[maxSize];		
	}

	public Euler2ndOrder(double[] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        n = (int)(Math.sqrt(migration_rates.length) + 0.5);        
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];tCR = new double[states]; sumDotStates = new double[states];   	
    	hasIndicators = false;
    	hasMultiplicator = false;
    	
    	iterations=0;
	}
	
	public Euler2ndOrder(double[] migration_rates, int[] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
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
    	sumStates = new double[states];tCR = new double[states];
    	hasIndicators = true;
    	hasMultiplicator = false;
    	
    	iterations=0;
	}

	public Euler2ndOrder(int[] multiplicator, double[] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        n = (int)(Math.sqrt(migration_rates.length) + 0.5);        
        this.multiplicator = multiplicator;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];tCR = new double[states]; sumDotStates = new double[states];   	
    	hasIndicators = false;
    	hasMultiplicator = true; 
    	
    	iterations=0;
	}
	
	public Euler2ndOrder(int[] multiplicator, double[] migration_rates, int[] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step) {
		this.max_step = max_step;
		this.epsilon = epsilon;
        this.migration_rates = migration_rates;
        this.multiplicator = multiplicator;
        n = (int)(Math.sqrt(migration_rates.length) + 0.5);        
        this.indicators = indicators;
        n2 = indicators.length / 2;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        this.dimension = this.lineages*this.states;
    	sumStates = new double[states];tCR = new double[states]; sumDotStates = new double[states];   	
    	hasIndicators = true;
    	hasMultiplicator = true;   
    	
    	iterations=0;
	}

	@Override
	public void calculateValues(double duration, double[] p, int length){
		double[] pDot = linProbs_tmpdt; 
		double[] pDotDot = linProbs_tmpddt; 
		double[] pDotDotDot = linProbs_tmpdddt;
		calculateValues(duration, p, pDot, pDotDot, pDotDotDot, length);
	}
	
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot){
		calculateValues(duration, p, pDot, pDotDot, pDotDotDot, pDot.length);
	}
	
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length){
		clearArray(pDotDot, length);
		clearArray(pDotDotDot, length);

		if (hasMultiplicator) {
			while (duration > 0){
		    	iterations++;
				//pDot = new double[length];
				clearArray(pDot, length);
				computeDerivativesWithMultiplicator(p, pDot, pDotDot, pDotDotDot, length);
				computeSecondDerivateWithMultiplicator(p, pDot, pDotDot, length);
				approximateThirdDerivate(pDotDot, pDotDotDot, length);
				duration = updateP(duration, p,  pDot, pDotDot, pDotDotDot, length - 1);
				
				if (iterations>10000){
					System.err.println("too many iterations, erturn negative infinity");
					p[length-1] = Double.NEGATIVE_INFINITY;
					break;
				}
			}			
		} else {
			while (duration > 0){
		    	iterations++;
				//pDot = new double[length];
				clearArray(pDot, length);
				computeDerivatives(p, pDot, pDotDot, pDotDotDot, length);
				computeSecondDerivate(p, pDot, pDotDot, length);
				approximateThirdDerivate(pDotDot, pDotDotDot, length);
				duration = updateP(duration, p,  pDot, pDotDot, pDotDotDot, length - 1);
				
				if (iterations>10000){
					System.err.println("too many iterations, erturn negative infinity");
					p[length-1] = Double.NEGATIVE_INFINITY;
					break;
				}
			}
		}			
	}	
	
	private void clearArray(double[] v, int n) {
		for (int i = 0; i < n; i++) {
			v[i] = 0.0;
		}		
	}

	private double updateP (double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length){
		final double max_dotdotdot = maxAbs(pDotDotDot, length);	
		
		//double timeStep = FastMath.min(FastMath.pow(epsilon*6/max_dotdotdot, C), FastMath.min(duration, max_step));

		double timeStep = FastMath.min(FastMath.cbrt(epsilon*6/max_dotdotdot), FastMath.min(duration, max_step));
		double timeStepSquare = timeStep * timeStep * 0.5;
		
		for (int i = 0; i < length; i++) {
			double new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
			double diff = FastMath.abs(new_val - p[i]);
			while (new_val > 1 || new_val < 0 || diff > 0.2) {
				timeStep *= 0.9;
				timeStepSquare = timeStep * timeStep * 0.5;
				new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
				diff = FastMath.abs(new_val - p[i]);
			}			
		}
		
		updateP2(timeStep, timeStepSquare, p, length + 1, pDot, pDotDot);
		
		// normalize to ensure stability
		for (int i = 0; i < lineages; i ++) {
			normalise(i, p);
		}
		
		duration -= timeStep;
		return duration;
	}
	
	
	private double maxAbs(double[] pDotDotDot, int length) {
		double max_dotdotdot = 0.0;
		for (int i = 0; i < length; i++) {
			max_dotdotdot = FastMath.max(max_dotdotdot, FastMath.abs(pDotDotDot[i]));
		}
		return max_dotdotdot;
	}

	    
	private void normalise(final int i, final double[] p) {
		final int k = states * i;
		double linSum = 0;
		
		int u = k;
		for (int j = 0; j < states; j++) {
			final double x = p[u++];
			linSum += x;
			if (x < 0.0) {
				System.err.println(Arrays.toString(p));
				System.exit(0);
			}
		}
		u = k;
		for (int j = 0; j < states; j++) {
			p[u++] /= linSum;
		}
	}

	private void updateP2(final double timeStep, final double timeStepSquare, final double[] p, final int length, final double[] pDot,
			final double[] pDotDot) {
		for (int i = 0; i < length; i++) {
			p[i] += pDot[i] * timeStep
			     + pDotDot[i] * timeStepSquare;
		}
	}

	public void computeDerivatives (double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
		
    	double migrates;
    	// Compute the sum of line state probabilities for each state
     	clearArray(sumStates, states);
     	calcSumStates(sumStates, p);
     	
    	// Calculate the change in the lineage state probabilities for every lineage in every state
		int currlin = 0;
    	for (int i = 0; i<lineages; i++){
    		
    		double sumCoal = 0;
    		int k = currlin;
    		for (int j = 0; j < states; j++) {
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
     		pDot[length-1] -= sumCoal;
     		
     		k = currlin;
    		for (int j = 0; j < states; j++) {    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
    		currlin += states;
    	}
    	
    	
    	// Calculate the probability of a lineage changing states
    	if (hasIndicators){
			for (int j = 0; j < indicators.length/2; j++){
				int source = indicators[j * n2 + 0];
				int sink = indicators[j * n2 + 1];
				double mrate = migration_rates[source * n + sink];
				int k = source;
				int m = sink;
		    	for (int i = 0; i<lineages; i++){
					migrates = p[k] * mrate;
					pDot[m] += migrates;
					pDot[k] -= migrates;  			
		    		k += states;
		    		m += states;
		    	}
	    	}    
    	} else {
    		int u = 0;
        	for (int i = 0; i < lineages; i++) {
        		// Calculate the probability of a lineage changing states
        		for (int j = 0; j < states; j++) {
            		int v = u;
        			double pj = p[u];
        			for (int k = j + 1; k < states; k++){    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			}// j
        			u++;
        		}// j
        	}// lineages       		
    	}
    	
		pDot[length-1]  /= 2;

    }
        
    private void calcSumStates(final double [] sumStates, final double[] p) {
     	int u = 0;
    	for (int i = 0; i < lineages; i++) {
    		for (int j = 0; j < states; j++) {
				sumStates[j] += p[u++];
    		}
    	}
	}
    
	public void computeSecondDerivate (double[] p, double[] pDot, double[] pDotDot, int length){
    	clearArray(sumDotStates, states);
    	calcSumStates(sumDotStates, pDot);
	
    	// Calculate the change in the lineage state probabilities for every lineage in every state
		int currlin = 0;
    	for (int i = 0; i < lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		int k = currlin;
    		for (int j = 0; j < states; j++) {
    			//final double pk = p[k];
    			//final double pDotk = pDot[k];
//    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - p[k]) + p[k] * (sumDotStates[j] - pDot[k]));
//    			pCoalRate += coalescent_rates[j] * (pDotk * (sumStates[j] - pk) + pk * (sumDotStates[j] - pDotk));
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
    		
//    		k = currlin;
//    		for (int j = 0; j<states; j++) {
//    			pDotDot[k] *= pDot[k];
//    			k++;
//    		}

    		k = currlin;
    		for (int j = 0; j < states; j++) {
//    			final double pDotk = pDot[k];
//    			pDotDot[k] = pDotDot[k] * pDotk + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDotk));
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j

    		pDotDot[length-1] -= pCoalRate;

    		currlin += states;
    	}// lineages 
    	
		double migrates;
		
		
		// Calculate the probability of a lineage changing states
		if (hasIndicators){
			for (int j = 0; j < indicators.length/2; j++){
				int source = indicators[j * n2 + 0];
				int sink = indicators[j * n2 + 1];
				double mrate = migration_rates[source * n + sink];
		    	for (int i = 0; i<lineages; i++){
					migrates = pDot[states*i+source]*mrate;
					pDotDot[states*i+sink] += migrates;
					pDotDot[states*i+source] -= migrates;  			
		    	}
	    	}    
		}else{
			int u = 0;
        	for (int i = 0; i<lineages; i++){
        		// Calculate the probability of a lineage changing states
        		for (int j = 0; j < states; j++){
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j    	
        			u++;
        		}// j
        	}// lineages    
			
		}
		pDotDot[length-1] /= 2;
    }
        
	public void approximateThirdDerivate (double[] pDotDot, double[] pDotDotDot, int length) {	
    	double migrates;
    	
    	// Calculate the change in the lineage state probabilities for every lineage in every state
		for (int u = 0; u < length - 1; u++) {
			pDotDotDot[u] *= pDotDot[u];
		}
    	
		// Calculate the probability of a lineage changing states
    	if (hasIndicators){
			for (int j = 0; j < indicators.length/2; j++){
				int source = indicators[j * n2 + 0];
				int sink = indicators[j * n2 + 1];
				double mrate = migration_rates[source * n + sink];
		    	for (int i = 0; i<lineages; i++){
					migrates = pDotDot[states * i + source] * mrate;
					pDotDotDot[states * i + sink] += migrates;
					pDotDotDot[states * i + source] -= migrates;  			
		    	}
	    	}    
    	} else {
			for (int j = 0; j < states; j++){
				for (int k = 0; k < states; k++){  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} 			
				}
	    	}
    	}
    }

	public void computeDerivativesWithMultiplicator(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
		
    	double migrates;
    	// Compute the sum of line state probabilities for each state
     	clearArray(sumStates, states);
     	{
	     	int u = 0;
	    	for (int i = 0; i<lineages; i++) {
	    		for (int j = 0; j<states; j++) {
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
	    	}
     	}
    		
    	// Calculate the change in the lineage state probabilities for every lineage in every state
		//double[] tCR =  new double[states];
    	for (int i = 0; i<lineages; i++){
    		double sumCoal = 0;
    		int currlin = states*i;
    		for (int j = 0; j<states; j++){
    			tCR[j] = coalescent_rates[j] *  (sumStates[j] - p[currlin+j]);
    			sumCoal += p[currlin+j]*tCR[j];
    		}
       		pDot[length-1] -= multiplicator[i]*sumCoal;
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
			for (int j = 0; j < indicators.length/2; j++){
				int source = indicators[j * n2 + 0];
				int sink = indicators[j * n2 + 1];
				double mrate = migration_rates[source * n + sink];
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
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j    			 
        		}// j
        	}// lineages       		
    	}
    	
		pDot[length-1]  /= 2;

    }
        
    public void computeSecondDerivateWithMultiplicator(double[] p, double[] pDot, double[] pDotDot, int length){  
    	//double[] sumDotStates = new double[states];;
    	clearArray(sumDotStates, states);
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];   
	
    	// Calculate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		int currlin = states*i;
    		for (int j = 0; j<states; j++)
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));     
    		
    		for (int j = 0; j<states; j++)
    			pDotDot[currlin+j] *= pDot[currlin+j];
    		
   			pDotDot[length-1] -= multiplicator[i]*pCoalRate;

    		for (int j = 0; j<states; j++){
    			pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    		}// j    		
    	}// lineages 
    	
		double migrates;
		
		
		// Calculate the probability of a lineage changing states
		if (hasIndicators){
			for (int j = 0; j < indicators.length/2; j++){
				int source = indicators[j * n2 + 0];
				int sink = indicators[j * n2 + 1];
				double mrate = migration_rates[source * n + sink];
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
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j    			 
        		}// j
        	}// lineages    
			
		}
		pDotDot[length-1] /= 2;
    }
    
 }


