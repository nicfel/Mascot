package beast.mascot.dynamics;


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;


@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class Skyline extends Dynamics  {

    public Input<RealParameter> NeInput = new Input<>("Ne", "input of effective population sizes", Validate.REQUIRED);    
    public Input<RealParameter> b_mInput = new Input<>("backwardsMigration", "input of backwards in time migration rates");    
    public Input<RealParameter> f_mInput = new Input<>("backwardsMigration", "input of backwards in time migration rates", Validate.XOR, b_mInput);    
    public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts", "time points at which rates are changed (Default no rate shifts)", Validate.REQUIRED);    
    public Input<Tree> treeInput = new Input<>("tree", "the timing of rate shifts is relative to the tree height (Default true)", Validate.OPTIONAL);

	private boolean isBackwardsMigration;
    

    @Override
    public void initAndValidate() {
    	if (b_mInput.get()!=null)
    		isBackwardsMigration = true;
    	else
    		isBackwardsMigration = false;
    }


    /**
     * Returns the time to the next interval.
     */
    public double getInterval(int i) {
    	if (treeInput.get()!=null)
    		return rateShiftsInput.get().getArrayValue(i) * treeInput.get().getRoot().getHeight();
    	else
    		return rateShiftsInput.get().getArrayValue(i);
    }   
    
    public boolean intervalIsDirty(int i){
    	boolean intervalIsDirty = false;   	
    	
    	// Check if rate shifts, Ne's or migration rates have changed for the current interval    	
    	if (rateShiftsInput.get() != null)
    		if (rateShiftsInput.get().isDirty(i))
    			intervalIsDirty = true;
    	
    	for (int j = i*dimensionInput.get(); j < (i+1)*dimensionInput.get(); j++)
    		if (NeInput.get().isDirty(j))
    			intervalIsDirty = true;  
    
    	if (isBackwardsMigration){
			for (int j = i*dimensionInput.get()*(dimensionInput.get()-1); j < (i+1)*dimensionInput.get()*(dimensionInput.get()-1); j++)
				if (b_mInput.get().isDirty(j))
					intervalIsDirty = true; 
    	}else{
			for (int j = 0; j < f_mInput.get().getDimension(); j++)
				if (f_mInput.get().isDirty(j))
					intervalIsDirty = true;    		
    	}		
    	return intervalIsDirty;

    }   
 
	@Override
    public double[] getCoalescentRate(int i){
    	double[] Ne = new double[dimensionInput.get()];
    	int c = 0;
    	for (int j = i*dimensionInput.get(); j < (i+1)*dimensionInput.get(); j++){
    		Ne[c] = 1/NeInput.get().getArrayValue(j);
    		c++;
    	}
    	
    	return Ne;

    }
    
	@Override    
    public double[][] getBackwardsMigration(int i){
    	double[][] m = new double[dimensionInput.get()][dimensionInput.get()];
    	
    	if (isBackwardsMigration){
	    	int c = i*dimensionInput.get()*(dimensionInput.get()-1);
	    	for (int a = 0; a < dimensionInput.get(); a++){
	    		for (int b = 0; b < dimensionInput.get(); b++){
	    			if (a!=b){
	    				m[a][b] = b_mInput.get().getArrayValue(c);
	    				c++;
	    			}
	    		}
	    	}
    	}else{
	    	int c = 0;
	    	for (int a = 0; a < dimensionInput.get(); a++){
	    		for (int b = 0; b < dimensionInput.get(); b++){
	    			if (a!=b){
	    				m[a][b] = f_mInput.get().getArrayValue(c)
	    						*NeInput.get().getArrayValue(i*dimensionInput.get() + b)
	    							/NeInput.get().getArrayValue(i*dimensionInput.get() + a);
	    				c++;
	    			}
	    		}
	    	}    		
    	}
    	
    	return m;  	
    }  	

	@Override
	public void recalculate() {
	}



    
}