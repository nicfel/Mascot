
package beast.mascot.dynamics;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

@Description("Probabilistic representation that can produce " +
        "a log probability for instance for running an MCMC chain.")
public abstract class Dynamics extends CalculationNode  {
	
    public Input<Integer> dimensionInput = new Input<>("dimension", "the number of different states ", Validate.REQUIRED);

    private boolean dynamicsKnown;

    /**
     * recalculate the dynamics
     */
    public abstract void recalculate();
    
    /**
     * get the time to the next rate shift
     */    
	public abstract double getInterval(int i);

    /**
     * check if next interval is dirty, i.e. parameters changed in that interval
     */    
	public abstract boolean intervalIsDirty(int i);
	
    /**
     * get the effective population sizes for the next interval
     */    
	public abstract double[] getCoalescentRate(int i);
	
    /**
     * get the migration rates for the next interval
     */    
	public abstract double[][] getBackwardsMigration(int i);


    @Override
    public void initAndValidate() {
        // nothing to do
    }
    
    /**
     * Intended to be overridden by stochastically estimated distributions.
     * Used to disable target distribution consistency checks implemented in
     * the MCMC class which do not apply to stochastic distributions.
     * 
     * @return true if stochastic.
     */
    public boolean isStochastic() {
        return false;
    }
    
    @Override
	protected boolean requiresRecalculation(){
    	dynamicsKnown = false; 	
    	return dynamicsKnown;
    }
    
    public boolean areDynamicsKnown(){
    	return dynamicsKnown;
    }

	public int getDimension(){
		return dimensionInput.get();
	}
	

} // class Distribution
