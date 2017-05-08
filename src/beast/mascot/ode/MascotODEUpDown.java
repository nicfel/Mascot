package beast.mascot.ode;

import java.util.Arrays;
import org.apache.commons.math4.ode.FirstOrderDifferentialEquations;

/**
 * @author Nicola Felix Mueller
 */
public class MascotODEUpDown implements FirstOrderDifferentialEquations {
    
    
	double[][] migration_rates;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
    
    boolean belowzero = false;

    // MascotODEUpDown
    public MascotODEUpDown(double[][] migration_rates, double[] coalescent_rates, int lineages , int states){
        this.migration_rates = migration_rates;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        belowzero = false;
        this.dimension = this.lineages*this.states;
    }
    
    
    public int getDimension() {
    	// Each lineage has additionally nr states entries that are used 
    	// to compute the transition probabilities
        return dimension + (states*states*lineages);
    }
    
    public void computeDerivatives(double t, double[] p, double[] pDot) {
    	// normalize each lineage such that the probability of the lineage not having coalesced
    	// is 1. This is done such that the coalescent rate is calculated conditional on all 
    	// lineages still existing
    	double[] p_norm = normalizeLineages(p);
  	
    	// Compute the sum of line state probabilities for each state
    	double[] sumStates = getSumStates(p_norm);
    	double dTdt = 0.0;
    	double[] dTdtStates = new double[states];   	
    	
    	for (int i = 0; i<lineages; i++){
    		for (int j = 0; j<states; j++){
    			dTdtStates[j] += p_norm[i*states+j]*(sumStates[j] - p_norm[i*states+j]);
    		}
    	}

		for (int j = 0; j<states; j++){
			dTdt += coalescent_rates[j]*dTdtStates[j];
		}

    	// Calculate pDot from migration
		pDot = getMigration(p, pDot); 	

    	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	pDot = getCoalescence(p, p_norm, pDot, sumStates, dTdt);
//    	System.out.println(Arrays.toString(pDot));
    	
		pDot = getMigrationTransition(p, pDot);
    	
    	// Calculate the transition probabilities due to coalescence
    	pDot = getCoalescenceTransition(p, p_norm, pDot, sumStates, dTdt);   
    }
    
    private double[] normalizeLineages(double[] p_in){
    	double[] p_norm = new double[p_in.length];
    	for (int i = 0; i<lineages; i++){
    		double linSum=0.0;
    		for (int j = 0; j<states; j++){
    			linSum += p_in[states*i+j];    			
    		}
    		for (int j = 0; j<states; j++){
    			p_norm[states*i+j] = p_in[states*i+j]/linSum;   
    		}
    	}
    	return p_norm;    	
    }
    
 
    private double[] getSumStates(double[] p_norm){
    	double[] sumStates = new double[states];
    	for (int j = 0; j<states; j++)
    		sumStates[j] = 0;
    	
    	int c = 0;    	
    	for (int i = 0; i<lineages; i++){
    		for (int j = 0; j<states; j++){
    			sumStates[j] += p_norm[c];    			
    			c++;
    		}
    	}
    	return sumStates;
    }
    private double[] getMigration(double[] p, double[] pDot){
//    	pDot = new double[p.length];
		for (int j = 0; j<states; j++)
			for (int i = 0; i<lineages; i++)
				pDot[states*i+j] = 0.0;
			
			
		for (int j = 0; j<states; j++){
			for (int k = 0; k<states; k++){
				if (k != j){
					for (int i = 0; i<lineages; i++){
						double migration = p[states*i+j]*migration_rates[j][k];
						pDot[states*i+j] -= migration;
						pDot[states*i+k] += migration;
					}//i
				}// if k!=j
			}// k
		}// j
    	return pDot;
    }
    private double[] getMigrationTransition(double[] p, double[] pDot){
    	for (int j = 0; j < states; j++)
    		for (int k = 0; k < states; k++)
    			for (int i = 0; i < lineages; i++)
					pDot[lineages*states + i*states*states + j*states + k] = 0.0;
    	for (int j = 0; j < states; j++){ // initial value
    		for (int k = 0; k < states; k++){ // which current state
    			for (int l = 0; l < states; l++){  // which flow to state
    				if (k != l){
	    				for (int i = 0; i < lineages; i++){ // which lineages
	    					double migration = p[lineages*states + i*states*states + j*states + k]
	    							* migration_rates[k][l];	    					
	    					
	    					pDot[lineages*states + i*states*states + j*states + k] -= migration;
							pDot[lineages*states + i*states*states + j*states + l] += migration;
	    				}// i
    				}// if
    			}// l
    		}// k
    	}// j
    	return pDot;
    }
  
    private double[] getCoalescence(double[] p, double[] p_norm, double[] pDot, double[] sumStates, double dTdt){
    	for (int i = 0; i<lineages; i++){
    		// get the change of P_t(T) due to coalescence w/o lineage i
    		double dTdtWoI = dTdt;
    		for (int j = 0; j < states; j++){
    			dTdtWoI -= 2*p_norm[i*states+j]*coalescent_rates[j]*(sumStates[j] - p_norm[i*states+j]);
    		}        	        	
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			
    			// Calculate the Derivate of p:
    			pDot[states*i+j] -=
    					  2*coalescent_rates[j] * p[states*i+j] * (sumStates[j] - p_norm[states*i+j]) -
    					  	p[states*i+j]*dTdtWoI;
    		}// j
    	}// lineages   
    	return pDot;
    }
    
    
    private double[] getCoalescenceTransition(double[] p, double[] p_norm, double[] pDot, double[] sumStates, double dTdt){
    	for (int i = 0; i<lineages; i++){
    		// get the change of P_t(T) due to coalescence w/o lineage i
    		double dTdtWoI = dTdt;
    		for (int j = 0; j < states; j++){
    			dTdtWoI -= 2*p_norm[i*states+j]*coalescent_rates[j]*(sumStates[j] - p_norm[i*states+j]);
    		}        	        	
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
        		for (int k = 0; k < states; k++){ // loop over the different states, only the active states   			
    			// Calculate the Derivate of p:
    			pDot[lineages*states + i*states*states + j*states + k] -=
    					  2*coalescent_rates[k] * p[lineages*states + i*states*states + j*states + k] * (sumStates[k] - p_norm[states*i+k]) -
    					  	p[lineages*states + i*states*states + j*states + k]*dTdtWoI;
        		}
    		}// j
    	}// lineages     	
    	
    	return pDot;
    }
}