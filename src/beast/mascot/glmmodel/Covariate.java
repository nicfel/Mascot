package beast.mascot.glmmodel;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import cern.colt.Arrays;


public class Covariate extends BEASTObject  {
    
	public Input<List<Double>> valuesInput = 
			new Input<>("value", "start value(s) for this parameter. If multiple values are specified, they should be separated by whitespace.", 
					new ArrayList<>(), beast.core.Input.Validate.REQUIRED, getMax().getClass());

	private Double[] values;
	private List<String> rawValues;
	private int dimension;
	public String isTimeDependent = "false";
	public Boolean transformed = false;
	
	public Covariate() {}
	
	public Covariate(Double[] values, String id) {
		this.values = new Double[values.length];
		System.arraycopy(values, 0, this.values, 0, values.length);
		this.ID = id;
	}
	
	public Covariate(List<String> rawValues, String id) {
		this.rawValues = new ArrayList<>(rawValues);
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
    
    // log standardizes values
    public void transform() {
    	if (transformed)
    		return;

    	Double[] log_vals = new Double[values.length];
    	for (int i = 0; i < values.length; i++) {
    		if (values[i]<=0.0)
    			throw new IllegalArgumentException("predictor " + getID() + " should be log transformed and standardized but contains values that are 0 or smaller");
    		log_vals[i] = Math.log(values[i]);
    	}
    	double mean=0;
    	for (int i = 0; i < log_vals.length; i++) 
    		mean+=log_vals[i];
    	
    	mean/=log_vals.length;

    	double std=0;
    	for (int i = 0; i < log_vals.length; i++)
    		std+= (log_vals[i]-mean)*(log_vals[i]-mean);
    	
    	std/=log_vals.length-1;
    	std = Math.sqrt(std);
    			
    	for (int i = 0; i < log_vals.length; i++) 
    		log_vals[i] = (log_vals[i]-mean)/std;
    	
		System.arraycopy(log_vals, 0, values, 0, log_vals.length);
		transformed = true;
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
    
    // from raw data to double values
    public void initMigrationFromRawValues(HashMap<String, Integer> traitToType, int nrIntervals) {
    	isTimeDependent = "false";
    	
    	// build int map for traits
    	int[][] indicesMap = new int[traitToType.size()][traitToType.size()];
    	int c=0;
    	for (int a=0;a<traitToType.size();a++) {
    		for (int b=0;b<traitToType.size();b++) {
    			if (a!=b) {
    				indicesMap[a][b] = c;
    				c++;
    			}    				
    		}  
		}
    	
    	// TODO, add possibility for time dependence
    	
    	
    	// get the header values
    	String[] headervals = rawValues.get(0).replace("\\s+", "").split(",");  
    	
    	if ((headervals.length-1)%traitToType.size()!=0) {
			isTimeDependent = "wrong number of entries ";
			return;
    	}
    	
    	int time_points = (headervals.length-1)/traitToType.size();
    	
    	if (time_points!=1 && time_points!=nrIntervals) {
			isTimeDependent = "wrong number of entries";
			return;
    	}
    	
    	int dimension = (traitToType.size()-1)*traitToType.size() * nrIntervals;
    	this.values = new Double[dimension]; 
    	
    	for (int a = 1; a < rawValues.size(); a++) {
    		String[] splitvals = rawValues.get(a).replace("\\s+", "").split(",");
    		for (int b = 1; b < splitvals.length; b++) {
    			// get the from and to values
    			int from = traitToType.get(splitvals[0]);
    			int to = traitToType.get(headervals[b].split(":")[0]);
    			
    			int timepoint = 0;
    			
    			// check if there is a time information
    			if (headervals[b].contains(":")) {
    				timepoint = Integer.parseInt(headervals[b].split(":")[1]);
    				isTimeDependent = "true";
    			}
    			
    			if (from!=to) {
	    			int index = indicesMap[from][to];
	    			this.values[index + timepoint*c] = Double.parseDouble(splitvals[b]);
    			}
    		}    		    		
    	}
    	
    	if (isTimeDependent.contains("false")) {
    		for (int timepoint = 1; timepoint < nrIntervals; timepoint++) {
    			int index=0;
    			int dims = traitToType.size()*(traitToType.size()-1);
    	    	for (int a = 0; a < traitToType.size(); a++) {
    	    		for (int b = 0; b < traitToType.size(); b++) {
    	    			if (a!=b) {
    	    				this.values[index + timepoint*dims] = this.values[index];
    	    				index++;
    	    			}    	    			
    	    		}
    	    	}
    		}
    	}
    	
    	
    	
    	// convert values to arraylist
    	List<Double> inputvals = new ArrayList<>();
    	for (int i = 0; i < this.values.length;i++)
    		inputvals.add(this.values[i]);
    	
    	// set the input values
    	valuesInput.set(inputvals);
    	
    }
        
    // from raw data to double values
    public void initNeFromRawValues(HashMap<String, Integer> traitToType, int nrIntervals) {
    	// check if the first line starts with a location, if not, the predictor is time dependent    	
    	if (traitToType.containsKey(rawValues.get(0).replace("\\s+", "").split(",")[0]))
        	isTimeDependent = "false";
    	else
    		isTimeDependent = "true";    	
    	
    	if (isTimeDependent.contentEquals("true")) {
        	// get the header values
        	String[] headervals = rawValues.get(0).replace("\\s+", "").split(",");  
        	int time_points = headervals.length-1;
        	this.values = new Double[traitToType.size() *time_points]; 
        	
        	if (time_points!=1 && time_points!=nrIntervals) {
        		isTimeDependent = "wrong number of entries";
				return;
        	}

        	
        	// skip first line
        	for (int a = 1; a < rawValues.size(); a++) {
        		String[] splitvals = rawValues.get(a).replace("\\s+", "").split(",");
        		for (int b = 1; b < splitvals.length; b++) {
        			// get the from and to values
        			int from = traitToType.get(splitvals[0]);
        			
        			int timepoint = -1;
        			
        			// check if there is a time information
    				timepoint = Integer.parseInt(headervals[b]);
        			
	    			this.values[from + timepoint*traitToType.size()] = Double.parseDouble(splitvals[b]);
        			
        		}    		    		
        	}
    	}else {
        	// skip first line
        	this.values = new Double[traitToType.size()*nrIntervals]; 

        	for (int a = 0; a < rawValues.size(); a++) {
        		String[] splitvals = rawValues.get(a).replace("\\s+", "").split(",");
        		for (int b = 1; b < splitvals.length; b++) {
        			// get the from and to values
        			int from = traitToType.get(splitvals[0]);       			        			
	    			this.values[from] = Double.parseDouble(splitvals[b]);        			
        		}    		    		
        	}
        	
    		for (int timepoint = 1; timepoint < nrIntervals; timepoint++) {
    	    	for (int a = 0; a < traitToType.size(); a++) {
    				this.values[a + timepoint*traitToType.size()] = this.values[a];
    	    	}
    		}       	
    	}    	

    	// convert values to arraylist
    	List<Double> inputvals = new ArrayList<>();
    	for (int i = 0; i < this.values.length;i++)
    		inputvals.add(this.values[i]);
    	
    	// set the input values
    	valuesInput.set(inputvals);    	
    }

}
