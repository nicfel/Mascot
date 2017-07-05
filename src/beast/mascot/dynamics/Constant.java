package beast.mascot.dynamics;


import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;


@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class Constant extends Dynamics  {

    public Input<RealParameter> NeInput = new Input<>("Ne", "input of effective population sizes", Validate.REQUIRED);    
    public Input<RealParameter> b_mInput = new Input<>("backwardsMigration", "input of backwards in time migration rates");    
    public Input<RealParameter> f_mInput = new Input<>("forwardsMigration", "input of backwards in time migration rates", Validate.XOR, b_mInput);    

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
    	return Double.POSITIVE_INFINITY;
    }   
    
    public boolean intervalIsDirty(int i){
    	boolean intervalIsDirty = false;  	    	
    	
    	for (int j = 0; j < dimensionInput.get(); j++)
    		if (NeInput.get().isDirty(j))
    			intervalIsDirty = true;  
    
    	if (isBackwardsMigration){
			for (int j = 0; j < b_mInput.get().getDimension(); j++)
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
    	double[] coal = new double[dimensionInput.get()];
    	int c = 0;
    	for (int j = 0; j < dimensionInput.get(); j++){
    		coal[c] = 1/(2*NeInput.get().getArrayValue(j));
    		c++;
    	}
    	return coal;

    }
    
	@Override    
    public double[][] getBackwardsMigration(int i){
    	double[][] m = new double[dimensionInput.get()][dimensionInput.get()];
    	
    	if (isBackwardsMigration){
	    	int c = 0;
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
	    						*NeInput.get().getArrayValue(b)
	    							/NeInput.get().getArrayValue(a);
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