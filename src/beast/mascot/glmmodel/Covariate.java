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
	public Boolean isTimeDependent = false;
	public Boolean transform = false;
	
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
    public void initMigrationFromRawValues(HashMap<String, Integer> traitToType) {
    	isTimeDependent = false;
    	
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
    		System.out.println(Arrays.toString(indicesMap[a]));
		}
    	
    	// TODO, add possibility for time dependence
    	
    	
    	// get the header values
    	String[] headervals = rawValues.get(0).replace("\\s+", "").split(",");  
    	
    	if ((headervals.length-1)%traitToType.size()!=0) {
    		System.err.println("incorrect dimension of predictors, values either missing, or too many of them");
    		return;
    	}
    	
    	int time_points = (headervals.length-1)/traitToType.size();
    	int dimension = (traitToType.size()-1)*traitToType.size() *time_points;
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
    				isTimeDependent = true;
    			}
    			
    			if (from!=to) {
	    			int index = indicesMap[from][to];
	    			this.values[index + timepoint*c] = Double.parseDouble(splitvals[b]);
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
    public void initNeFromRawValues(HashMap<String, Integer> traitToType) {
    	// check if the first line starts with a location, if not, the predictor is time dependent    	
    	if (traitToType.containsKey(rawValues.get(0).replace("\\s+", "").split(",")[0]))
        	isTimeDependent = false;
    	else
    		isTimeDependent = true;    	
    	
    	if (isTimeDependent) {
        	// get the header values
        	String[] headervals = rawValues.get(0).replace("\\s+", "").split(",");  
        	int time_points = headervals.length-1;
        	this.values = new Double[traitToType.size() *time_points]; 
        	
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
        	for (int a = 0; a < rawValues.size(); a++) {
        		String[] splitvals = rawValues.get(a).replace("\\s+", "").split(",");
        		for (int b = 1; b < splitvals.length; b++) {
        			// get the from and to values
        			int from = traitToType.get(splitvals[0]);
        			        			
	    			this.values[from ] = Double.parseDouble(splitvals[b]);        			
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
