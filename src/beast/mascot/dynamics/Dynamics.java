
package beast.mascot.dynamics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.TraitSet;

@Description("Probabilistic representation that can produce " +
        "a log probability for instance for running an MCMC chain.")
public abstract class Dynamics extends CalculationNode  {
	
    public Input<Integer> dimensionInput = new Input<>("dimension", "the number of different states." + 
    " if -1, it will use the number of different types ", -1);
    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set.  Used only by BEAUti.");


    public boolean hasIndicators = false;
    
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

	
	/**
	 * get indicator variables for migration
	 * @param i
	 * @return
	 */
	public int[][] getIndicators(int i){
		return null;
	}


	HashMap<String, Integer> traitToType = new HashMap<>(); 
	HashMap<Integer, String> reverseTraitToType; 

    @Override
    public void initAndValidate() {    	
    	if (typeTraitInput.get()!=null){
    		traitToType = new HashMap<>();
    		reverseTraitToType = new HashMap<>();
    		List<String> taxa = typeTraitInput.get().taxaInput.get().asStringList();
    		ArrayList<String> unique = new ArrayList<>();
    		for (int i = 0; i < taxa.size(); i++)
    			unique.add(typeTraitInput.get().getStringValue(taxa.get(i)));
    		
    		Collections.sort(unique);
    		for (int i = unique.size()-2; i > -1; i--)
    			if(unique.get(i+1).equals(unique.get(i)))
    				unique.remove(i+1);
    		  
    		for (int i = 0; i < unique.size(); i++)
    			traitToType.put(unique.get(i), i);
    		for (int i = 0; i < unique.size(); i++)
    			reverseTraitToType.put(i, unique.get(i));
    	}
    	
    	if (dimensionInput.get()>1 && dimensionInput.get()<traitToType.size())
            throw new IllegalArgumentException("dimension is not -1 (undefined) and smaller " +
            		"than the number of different traits");

    	
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
	

	public int getNrTypes(){
		return traitToType.size();
	}

	public int getValue(String id) {
		return traitToType.get(typeTraitInput.get().getStringValue(id));	
	}
	
	public String getStringStateValue(int state){
		return reverseTraitToType.get(state);
	}
	
	
} // class Distribution