package beast.mascot.ode;

import java.util.Arrays;

import org.apache.commons.math4.exception.MaxCountExceededException;
import org.apache.commons.math4.ode.FirstOrderDifferentialEquations;

import beast.core.Description;


/**
 * @author Nicola Felix Mueller
 */
@Description("Describes the derivates of lineage state probabilities for LISCO as described in Mueller et al., 2016")
public class MascotConditionalODE implements FirstOrderDifferentialEquations {

	double[][] migration_rates;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
    
    boolean belowzero = false;

    // constructor
    public MascotConditionalODE(double[][] migration_rates, double[] coalescent_rates, int lineages , int states){
        this.migration_rates = migration_rates;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        belowzero = false;
        this.dimension = this.lineages*this.states;
    }

    public int getDimension() {
        return dimension+1;
    }

    public void computeDerivatives (double t, double[] p, double[] pDot) {
    	// normalize each lineage such that the probability of the lineage not having coalesced
    	// is 1. This is done such that the coalescent rate is calculated conditional on all 
    	// lineages still existing/ not having coalesced
//    	double[] p_norm = normalizeLineages(p);
    	
    	// Compute the sum of line state probabilities for each state
    	double[] sumStates = new double[states];
     	
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


}