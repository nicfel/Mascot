package beast.mascot.dynamics;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;


@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class GLM extends Dynamics  {	
    
    public Input<List<String>> migrationCovariatesInput = new Input<>("migrationCovariates", "input of migration covariates", Validate.REQUIRED);
    public Input<RealParameter> migrationScalerInput = new Input<>("migrationScaler", "input of migration covariates scaler", Validate.REQUIRED);    
    public Input<BooleanParameter> migrationIndicatorInput = new Input<>("migrationIndicator", "input of migration covariates scaler");     
    
    public Input<List<String>> neCovariatesInput = new Input<>("NeCovariates", "input of migration covariates", Validate.REQUIRED);
    public Input<RealParameter> neScalerInput = new Input<>("migrationScaler", "input of migration covariates scaler", Validate.REQUIRED);    
    public Input<BooleanParameter> neIndicatorInput = new Input<>("migrationIndicator", "input of migration covariates scaler");    

    public Input<RealParameter> timeInput = new Input<>("time", "input of times at which to switch rates");    

    double[] times;
    
    //covar<time<dims>>
    List<List<Double[][]>> migrationCovariates;
    List<List<Double[]>> neCovariates;
    
    
    
    @Override
    public void initAndValidate() {
    	super.initAndValidate();
    	
    	// fill the time vector
    	times = new double[timeInput.get().getDimension()+1];
    	times[times.length-1] = Double.POSITIVE_INFINITY;
    	
    	for (int i = 0; i < (times.length-1); i++)
    		times[i] = timeInput.get().getArrayValue(i);
    	
    	List<Double[][]> addListMigration = new ArrayList<>();;
    	// read the migration Covatiates
    	migrationCovariates = new ArrayList<>();
    	for (int i = 0; i < migrationCovariatesInput.get().size(); i++){
    		addListMigration = new ArrayList<>();
    		if (!migrationCovariatesInput.get().get(i).contains(";")){
    			for (int j = 0; j < times.length; j++){
    				
    			}
   			
    		}else{
    			String[] tmp_string1 = migrationCovariatesInput.get().get(i).split(";");
    			for (int j = 0; j < tmp_string1.length; j++){
    				String[] tmp_string2 = tmp_string1[j].split(",");
    				Double[][] addCov = new Double[dimensionInput.get()][dimensionInput.get()];
    				int c = 0;
    				for (int a = 0; a < dimensionInput.get(); a++){
    					for (int b = 0; b < dimensionInput.get(); b++){
    						if (a != b){
    							addCov[a][b] = Double.parseDouble(tmp_string2[c].trim());
    						}//if
    					}//b
    				}//a
    				addListMigration.add(addCov);
    			}//j
    		}//if   
    		migrationCovariates.add(addListMigration);
    	}//i
		migrationCovariates.add(addListMigration);
		
		
    	List<Double[]> addListNe = new ArrayList<>();
    	// read the migration Covatiates
    	neCovariates = new ArrayList<>();
    	for (int i = 0; i < neCovariatesInput.get().size(); i++){
    		addListNe = new ArrayList<>();
    		if (!neCovariatesInput.get().get(i).contains(";")){
    			for (int j = 0; j < times.length; j++){
    				
    			}
   			
    		}else{
    			String[] tmp_string1 = neCovariatesInput.get().get(i).split(";");
    			for (int j = 0; j < tmp_string1.length; j++){
    				String[] tmp_string2 = tmp_string1[j].split(",");
    				Double[] addCov = new Double[dimensionInput.get()];
    				int c = 0;
    				for (int a = 0; a < dimensionInput.get(); a++){
						addCov[a] = Double.parseDouble(tmp_string2[c].trim());
    				}//a
    				addListNe.add(addCov);
    			}//j
    		}//if   
    		neCovariates.add(addListNe);
    	}//i
    	neCovariates.add(addListNe);
 	
    	
    }

    /**
     * Returns the time to the next interval.
     */
    public double getInterval(int i) {
    	return times[i];
    }   

    public boolean intervalIsDirty(int i){
    	boolean intervalIsDirty = false;  	    	
    	
    	for (int j = 0; j < migrationScalerInput.get().getDimension(); j++)
			if (migrationScalerInput.get().isDirty(j))
				intervalIsDirty = true; 
		
    	for (int j = 0; j < migrationIndicatorInput.get().getDimension(); j++)
			if (migrationIndicatorInput.get().isDirty(j))
				intervalIsDirty = true; 
		
    	for (int j = 0; j < neScalerInput.get().getDimension(); j++)
			if (neScalerInput.get().isDirty(j))
				intervalIsDirty = true; 
		
    	for (int j = 0; j < neIndicatorInput.get().getDimension(); j++)
			if (neIndicatorInput.get().isDirty(j))
				intervalIsDirty = true; 
	
    	return intervalIsDirty;

    }  
    
	@Override
    public double[] getCoalescentRate(int i){
    	double[] coal = new double[dimensionInput.get()];
    	for (int j = 0; j < dimensionInput.get(); j++){
    		coal[j] = 0;
    		for (int c = 0; c < neCovariates.size(); c++){
    			if (neIndicatorInput.get().getArrayValue(c) > 0.0)
    				coal[j] += neScalerInput.get().getArrayValue(c)*neCovariates.get(c).get(i)[j];
    		}
    	}
    	return coal;
    }
    
	@Override    
    public double[][] getBackwardsMigration(int i){
    	double[][] m = new double[dimensionInput.get()][dimensionInput.get()];
    	for (int a = 0; a < dimensionInput.get(); a++){
    		for (int b = 0; b < dimensionInput.get(); b++){
    			if (a != b){
	    			m[a][b] = 0;
		    		for (int c = 0; c < migrationCovariates.size(); c++){
		    			if (neIndicatorInput.get().getArrayValue(c) > 0.0)
		    				m[a][b] += migrationScalerInput.get().getArrayValue(c)*migrationCovariates.get(c).get(i)[a][b];
		    		}
    			}
    		}
    	}
    	return m;
    }

	@Override
	public void recalculate() {
		// TODO Auto-generated method stub
		
	}   

	

    
}