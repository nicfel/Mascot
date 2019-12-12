package beast.mascot.glmmodel;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import cern.colt.Arrays;


public class Covariate extends BEASTObject  {
    
	final public Input<List<Double>> valuesInput = 
			new Input<>("value", "start value(s) for this parameter. If multiple values are specified, they should be separated by whitespace.", 
					new ArrayList<>(), beast.core.Input.Validate.REQUIRED, getMax().getClass());

	final public Input<Boolean> logTransformInput = 
			new Input<>("logTransform", "if true, the values will be log transformed and standardized", false);

	final public Input<Boolean> timeInvariantInput = 
			new Input<>("timeInvariant", "if true, the predictor is assumed to be time invariant", false);
	
	private Double[] values;
	private int dimension;
	
	public Covariate() {}
	
	public Covariate(Double[] values, String id) {
		this.values = new Double[values.length];
		System.arraycopy(values, 0, this.values, 0, values.length);
		this.ID = id;
	}
	
	public Covariate(String id) {
		this.values = null;
		this.ID = id;
	}

	
    public void initAndValidate() {
    	Double[] valuesString = valuesInput.get().toArray((Double[]) Array.newInstance(getMax().getClass(), 0));
    	dimension = valuesString.length;
        values = (Double[]) Array.newInstance(getMax().getClass(), dimension);
        for (int i = 0; i < values.length; i++) {
            values[i] = valuesString[i % valuesString.length];
        }
    }

	
    public Double getValue() {
        return values[0];
    }

    public double getArrayValue(final int index) {
        return values[index];
    }

	public int getDimension() {
		return values.length;
	}
    
    Double getMax() {
        return Double.POSITIVE_INFINITY;
    }
    
    public String toString() {
    	return getID();
    }
    
}
