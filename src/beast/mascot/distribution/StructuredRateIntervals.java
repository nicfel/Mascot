package beast.mascot.distribution;


import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;


@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class StructuredRateIntervals extends CalculationNode  {

	// Inputs
    public Input<RealParameter> NeInput = new Input<>("Ne", "input of effective population sizes", Input.Validate.REQUIRED);    
    public Input<RealParameter> mInput = new Input<>("m", "input of backwards in time migration rates", Input.Validate.REQUIRED);    
    public Input<Integer> dimensionInput = new Input<>("dimension", "the number of different states ", Input.Validate.REQUIRED);
    public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts", "time points at which rates are changed (Default no rate shifts)");    
    public Input<Boolean> relativeRateShifts = new Input<>("relativeRateShifts", "the timing of rate shifts is relative to the tree height (Default true)", true);

	
	
    protected boolean intervalsKnown;   
    protected double intervalScaler = 1.0;
    

    @Override
    public void initAndValidate() {
    }

    /**
     * CalculationNode methods *
     */
    @Override
    protected boolean requiresRecalculation() {
        intervalsKnown = false;
        return true;
    }

    /**
     * Returns the time to the next interval.
     */
    public double getInterval(int i) {
    	if (rateShiftsInput.get()==null){
    		return Double.POSITIVE_INFINITY;
    	}else{
	    	if (relativeRateShifts.get())
	    		return rateShiftsInput.get().getArrayValue(i) * intervalScaler;
	    	else
	    		return rateShiftsInput.get().getArrayValue(i);
    	}
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
    
		for (int j = i*dimensionInput.get()*(dimensionInput.get()-1); j < (i+1)*dimensionInput.get()*(dimensionInput.get()-1); j++)
			if (mInput.get().isDirty(j))
				intervalIsDirty = true; 
		
    	return intervalIsDirty;
    }
   
    public void setIntervalScaler(double scaleFactor){
    	intervalScaler = scaleFactor;
    }
    
    public double[] getIntervalCoalRate(int i){
    	double[] Ne = new double[dimensionInput.get()];
    	int c = 0;
    	for (int j = i*dimensionInput.get(); j < (i+1)*dimensionInput.get(); j++){
    		Ne[c] = 1/NeInput.get().getArrayValue(j);
    		c++;
    	}
    	
    	return Ne;
    }
    
    public double[][] getIntervalMigRate(int i){
    	double[][] m = new double[dimensionInput.get()][dimensionInput.get()];
    	int c = i*dimensionInput.get()*(dimensionInput.get()-1);
    	for (int a = 0; a < dimensionInput.get(); a++){
    		for (int b = 0; b < dimensionInput.get(); b++){
    			if (a!=b){
    				m[a][b] = mInput.get().getArrayValue(c);
    				c++;
    			}
    		}
    	}
    	
    	return m;  	
    }
    
	public int getDimension() {
		return dimensionInput.get();
	}

    
    @Override
    protected void restore() {
        super.restore();        
    }

    @Override
    protected void store() {
        super.store();
    }


    
}