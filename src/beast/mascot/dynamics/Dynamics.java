
package beast.mascot.dynamics;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.TraitSet;

@Description("Probabilistic representation that can produce " +
        "a log probability for instance for running an MCMC chain.")
public abstract class Dynamics extends CalculationNode  {
	
    public Input<Integer> dimensionInput = new Input<>("dimension", "the number of different states ", 1);
    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set.  Used only by BEAUti.");


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

	TraitSet typeTraitSet; 

    @Override
    public void initAndValidate() {
    	
    	
//        if (typeTraitInput.get() != null)
//            typeTraitSet = typeTraitInput.get();

        // Construct type list.
        if (typeTraitSet == null) {
//            if (getTaxonset() != null) {
                TraitSet dummyTraitSet = new TraitSet();

                StringBuilder sb = new StringBuilder();
                for (int i=0; i<10; i++) {
                    if (i>0)
                        sb.append(",\n");
                    sb.append("lala").append("=NOT_SET");
                }
//                try {
//                    dummyTraitSet.initByName(
//                        "traitname", "type",
//                        "taxa", getTaxonset(),
//                        "value", sb.toString());
//                    dummyTraitSet.setID("typeTraitSet.t:"
//                        + BeautiDoc.parsePartition(getID()));
//                    setTypeTrait(dummyTraitSet);
//                } catch (Exception ex) {
//                    System.out.println("Error setting default type trait.");
//                }
//            }
        }

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
    
    public void setDynamicsKnown(){
    	dynamicsKnown = true;
    }

	public int getDimension(){
		return dimensionInput.get();
	}

	public int getValue(String id) {
		return (int) typeTraitInput.get().getValue(id);
			
	}
	
	
} // class Distribution
