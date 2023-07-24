package mascot.dynamics;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;

public class RateShifts extends CalculationNode implements Loggable {

    final public Input<String> dateTimeFormatInput = new Input<>("dateFormat", "the date/time format to be parsed, (e.g., 'dd/M/yyyy')", "dd/M/yyyy");
	public Input<List<Double>> valuesInput = new Input<>(
    		"value", "input of timings of rate shifts relative to the most recent sample", new ArrayList<>(), beast.base.core.Input.Validate.REQUIRED, getInputClass());  
    public Input<Tree> treeInput = new Input<>(
    		"tree", "input of tree to calculate the timings of the rates shifts relative to the tree height");       
    	
    boolean isRelative;
    public Double[] rateShifts;
    Tree tree;
    
    // first element which is larger than 0
	int timeIndexOffset;
	
	public RateShifts() {}

    
	@Override
	public void initAndValidate() {
		if (treeInput.get()==null) isRelative = false;		
		else isRelative = true;
		
		if (isRelative){ 
			tree = treeInput.get();
		}	
   		rateShifts = valuesInput.get().toArray((Double[]) Array.newInstance(getInputClass(), 0));
	}
	

	public double getValue(int i){
    	if (i >= rateShifts.length){
     		return Double.POSITIVE_INFINITY;
     	}

		if (isRelative) {
			return tree.getRoot().getHeight()*rateShifts[i];
		}else{
			return rateShifts[i];
		}
	}
	
	// get's the midpoint of the current intervals
	public double getIntervalMidpoint(int i){
    	if (i >= rateShifts.length-timeIndexOffset){
     		return Double.POSITIVE_INFINITY;
     	}

		if (isRelative) {
			double treeHeight = tree.getRoot().getHeight();
			if (i==0)
				return treeHeight*rateShifts[i]/2;
			else
				return treeHeight*(rateShifts[i]+rateShifts[i-1])/2;
		}else{
			if (i==0)
				return rateShifts[i]/2;
			else	
				return (rateShifts[i]+rateShifts[i-1])/2;
					
		}
	}
	
	// get's in which interval the current time is
	public int getIntervalNumber(double t){
		if (isRelative) {
			double treeHeight = tree.getRoot().getHeight();
			for (int i = 0; i < rateShifts.length; i++)
				if (t < treeHeight*rateShifts[i])
					return i;
		}else{
			for (int i = 0; i < rateShifts.length; i++)
				if (t < rateShifts[i])
					return i;					
		}
		return -1;
	}


	
	public int getCorrectedDimension(){
		return rateShifts.length;
	}
	
	public int getDimension(){
		if (rateShifts==null)
			return 1;
		
		return rateShifts.length;
	}
	
	public boolean somethingIsDirty(){		
		if (treeInput.get()!=null)
			return treeInput.get().somethingIsDirty();
		else
			return false;
	}

    private Class<?>  getInputClass() {
   		return ((Double) 0.0).getClass();
    }


	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < rateShifts.length; i++) {
			out.print("RateShift."+i + "\t");
			
		}
		
	}


	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < rateShifts.length; i++) {
			double val = 0;
			if (isRelative) {
				val = tree.getRoot().getHeight()*rateShifts[i];
			}else{
				val= rateShifts[i];
			}

			out.print(val + "\t");
			
		}

		
	}


	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
    
	

}

