package mascot.ode;


import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;

public class Euler2ndOrder27 extends Euler2ndOrder {
	double updateP (double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length){
		final double max_dotdotdot = maxAbs(pDotDotDot, length);	
		
		//double timeStep = FastMath.min(FastMath.pow(epsilon*6/max_dotdotdot, C), FastMath.min(duration, max_step));

		double timeStep = FastMath.min(FastMath.cbrt(epsilon*6/max_dotdotdot), FastMath.min(duration, max_step));
		double timeStepSquare = timeStep * timeStep * 0.5;
		
		for (int i = 0; i < length; i++) {
			double new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
			double diff = FastMath.abs(new_val - p[i]);
			int its = 0;
			while (new_val > 1 || new_val < 0 || diff > 0.2) {
				timeStep *= 0.9;
				timeStepSquare = timeStep * timeStep * 0.5;
				new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
				diff = FastMath.abs(new_val - p[i]);
				its++;

				if (its > 10000) {
//					System.err.println("cannot find proper time step, skip these parameter values");
					p[length-1] = Double.NEGATIVE_INFINITY;
					break;					
				}
			}			
		}
		
		if (p[length-1]==Double.NEGATIVE_INFINITY)
			return 0.0;

		updateP2(timeStep, timeStepSquare, p, length + 1, pDot, pDotDot);
		
		// normalize to ensure stability
		for (int i = 0; i < lineages; i ++) {
			normalise(i, p);
		}
		
		duration -= timeStep;
		return duration;
	}
	
	
	double maxAbs(double[] pDotDotDot, int length) {
		double max_dotdotdot = 0.0;
		for (int i = 0; i < length; i++) {
			max_dotdotdot = FastMath.max(max_dotdotdot, FastMath.abs(pDotDotDot[i]));
		}
		return max_dotdotdot;
	}

	    
	void normalise(final int i, final double[] p) {
		final int k = states * i;
		double linSum = 0;
		
		int u = k;
		
		double x;
		
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
{
			x = p[u++];
			linSum += x;
			if (x < 0.0) {
				p[p.length-1] = Double.NEGATIVE_INFINITY; 
				return;
			} // XXX
			
		}
		u = k;
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
{
			p[u++] /= linSum;
		}
	}

	void bailout(double[] p) {
		System.err.println(Arrays.toString(p));
		System.exit(0);		
	}
	
	void updateP2(final double timeStep, final double timeStepSquare, final double[] p, final int length, final double[] pDot,
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
		int currlin = 0, j, k;
    	for (int i = 0; i<lineages; i++){
    		
    		double sumCoal = 0;
    		k = currlin;
j = 0;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
j++;
{
    			tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
    			sumCoal += p[k] * tCR[j];
    			k++;
    		}
     		pDot[length-1] -= sumCoal;
     		
     		k = currlin;
j = 0;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[k] = coal;
    			pDotDotDot[k] = coal;
    			pDot[k] +=	p[k] * coal;
    			k++;
    		} // j
j++;
{    			
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
			for (j = 0; j < indicators.length/2; j++){
				int source = indicators[j * n2 + 0];
				int sink = indicators[j * n2 + 1];
				double mrate = migration_rates[source * n + sink];
				k = source;
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
    		int u = 0 , v; double pj;
        	for (int i = 0; i < lineages; i++) {
        		// Calculate the probability of a lineage changing states
j = 0;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
j++;
{
            		v = u;
        			pj = p[u];
        			for (k = j + 1; k < states; k++) {    
    					v++;
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[v] * migration_rates[k * n + j] -
    							   pj   * migration_rates[j * n + k];
    					pDot[u] += migrates;
    					pDot[v] -= migrates;
        			} // j XXX
        			u++;
        		}// j
        	}// lineages       		
    	}
    	
		pDot[length-1]  /= 2;

    }
        
    private void calcSumStates(final double [] sumStates, final double[] p) {
     	int u = 0, j;
    	for (int i = 0; i < lineages; i++) {
j = 0;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
j++;
{
				sumStates[j] += p[u++];
    		}
    	}
	}
    
	public void computeSecondDerivate (double[] p, double[] pDot, double[] pDotDot, int length){
    	clearArray(sumDotStates, states);
    	calcSumStates(sumDotStates, pDot);
	
    	// Calculate the change in the lineage state probabilities for every lineage in every state
		int currlin = 0, j;
    	for (int i = 0; i < lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		int k = currlin;
j = 0;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
    			k++;
    		}
    		

    		k = currlin;
j = 0;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j
j++;
{
    			pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
    			k++;
    		}// j

    		pDotDot[length-1] -= pCoalRate;

    		currlin += states;
    	}// lineages 
    	
		double migrates;
		
		
		// Calculate the probability of a lineage changing states
		if (hasIndicators){
			for (j = 0; j < indicators.length/2; j++){
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
j = 0;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
        			u++;
        		}// j
j++;
{
        			double pj = pDot[u];
        			int v = u + 1;
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[v]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[u] += migrates;
    					pDotDot[v] -= migrates;
    					v++;
        			}// j XXX
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
    		int k;
			for (int j = 0; j < states; j++){
k = 0;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
k++;
{  
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
			    	for (int i = 0; i<lineages; i++){
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
	    			} // XXX
			    	
				}
	    	}
    	}
    }

	public void computeDerivativesWithMultiplicator(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
		int j;
		
    	double migrates;
    	// Compute the sum of line state probabilities for each state
     	clearArray(sumStates, states);
     	{
	     	int u = 0;
	    	for (int i = 0; i<lineages; i++) {
j = 0;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
j++;
{
					sumStates[j] += multiplicator[i]*p[u++]; 
	    		}
	    	}
     	}
    		
    	// Calculate the change in the lineage state probabilities for every lineage in every state
		//double[] tCR =  new double[states];
    	for (int i = 0; i<lineages; i++){
    		double sumCoal = 0;
    		int currlin = states*i;
    		for (j = 0; j < states; j++){
    			tCR[j] = coalescent_rates[j] *  (sumStates[j] - p[currlin+j]);
    			sumCoal += p[currlin+j]*tCR[j];
    		}
       		pDot[length-1] -= multiplicator[i]*sumCoal;
    		for (j = 0; j < states; j++){    			
    			// Calculate the Derivate of p:
    			double coal = sumCoal - tCR[j];
    			pDotDot[currlin+j] = coal;
    			pDotDotDot[currlin+j] = coal;
    			pDot[currlin+j] +=	p[currlin+j] * coal;
    		}// j

    	}
    	
    	
    	// Calculate the probability of a lineage changing states
    	if (hasIndicators){
			for (j = 0; j < indicators.length/2; j++){
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
    		double pj;
        	for (int i = 0; i<lineages; i++){
        		int currlin = states*i;
        		// Calculate the probability of a lineage changing states
j = 0;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
j++;
{
        			pj = p[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = p[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDot[currlin+j] += migrates;
    					pDot[currlin+k] -= migrates;
        			}// j XXX
        			
        		}// j
        	}// lineages       		
    	}
    	
		pDot[length-1]  /= 2;

    }
        
    public void computeSecondDerivateWithMultiplicator(double[] p, double[] pDot, double[] pDotDot, int length){
    	int j;
    	//double[] sumDotStates = new double[states];;
    	clearArray(sumDotStates, states);
    	for (int i = 0; i<lineages; i++) {
j = 0;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
j++;
{
    			sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
    		}
    	}
	
    	// Calculate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){    		
    		double pCoalRate = 0.0;    		
    		int currlin = states*i;
j = 0;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
j++;
{
    			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));
    		}
    		
    		for (j = 0; j<states; j++) {
    			pDotDot[currlin+j] *= pDot[currlin+j];
    		}
    		
   			pDotDot[length-1] -= multiplicator[i]*pCoalRate;

    		for (j = 0; j<states; j++) {
    			pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
    		}// j    		
    	}// lineages 
    	
		double migrates;
		
		
		// Calculate the probability of a lineage changing states
		if (hasIndicators){
			for (j = 0; j < indicators.length/2; j++){
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
			double pj;
        	for (int i = 0; i<lineages; i++){
        		int currlin = states*i;
        		// Calculate the probability of a lineage changing states
j = 0;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
j++;
{
        			pj = pDot[currlin+j];
        			for (int k = j+1; k < states; k++){    
        				
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates = pDot[currlin+k]*migration_rates[k * n + j] -
    							pj*migration_rates[j * n + k];
    					pDotDot[currlin+j] += migrates;
    					pDotDot[currlin+k] -= migrates;
        			}// j  XXX
        			
        		}// j
        	}// lineages    
			
		}
		pDotDot[length-1] /= 2;
    }
    
 }


