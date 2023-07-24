package mascot.ode;


import beast.base.core.Description;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


/**
 * @author Nicola Felix Mueller
 */
@Description("Describes the derivates of lineage state probabilities for LISCO as described in Mueller et al., 2016")
public class MascotODE implements FirstOrderDifferentialEquations {

	double[] migration_rates;
	int n; // dimension of migration rate matrix
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
    
    boolean belowzero = false;

    // constructor
    public MascotODE(double[] migration_rates, double[] coalescent_rates, int lineages , int states){
        this.migration_rates = migration_rates;
        n = (int) (Math.sqrt(migration_rates.length) + 0.5);
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        belowzero = false;
        this.dimension = this.lineages*this.states;
    }

    public int getDimension() {
        return dimension;
    }

    public void computeDerivatives(double t, double[] p, double[] pDot) {
    	// normalize each lineage such that the probability of the lineage not having coalesced
    	// is 1. This is done such that the coalescent rate is calculated conditional on all 
    	// lineages still existing/ not having coalesced
    	double[] p_norm = normalizeLineages(p);
    	
    	// Compute the sum of line state probabilities for each state
    	double[] sumStates = new double[states];
    	
    	for (int i = 0; i<lineages; i++)
    		for (int j = 0; j<states; j++)
    			sumStates[j] += p_norm[states*i+j];    			
    		
    	// Calculate dTdt
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
		
//		pDot = new double[p.length];
    	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){
    		// get the change of P_t(T) due to coalescence w/o lineage i
    		double dTdtWoI = dTdt;
    		for (int j = 0; j < states; j++){
    			dTdtWoI -= 2*p_norm[i*states+j]*coalescent_rates[j]*(sumStates[j] - p_norm[i*states+j]);
    		}        	        	
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			double migrates = 0.0;
    			for (int k = 0; k<states; k++){
    				if (j != k){
    					    					
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates += p[states*i+k]*migration_rates[k * n + j] -
    									p[states*i+j]*migration_rates[j * n + k];
    				}
    			}// j    			 
    			
    			// Calculate the Derivate of p:
    			pDot[states*i+j] = migrates -
    					  2*coalescent_rates[j] * p[states*i+j] * (sumStates[j] - p_norm[states*i+j]) -
    					  	p[states*i+j]*dTdtWoI;
    		}// j
    	}// lineages   
    	
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


}