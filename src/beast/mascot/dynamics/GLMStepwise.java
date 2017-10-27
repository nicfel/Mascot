package beast.mascot.dynamics;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import com.google.common.primitives.Chars;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.mascot.glmmodel.GLMStepwiseModel;


@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class GLMStepwise extends Dynamics implements Loggable {	
    
	public Input<GLMStepwiseModel> migrationGLMInput = new Input<>(
			"migrationGLM", "input of migration GLM model", Validate.REQUIRED);
    
	public Input<GLMStepwiseModel> NeGLMInput = new Input<>(
			"NeGLM", "input of migration GLM model", Validate.REQUIRED);
    
    public Input<RealParameter> rateShiftsInput = new Input<>(
    		"rateShifts", "input of timings of rate shifts relative to the most recent sample", Validate.REQUIRED);    

	public Input<String> typesInput = new Input<>(
			"types", "input of the different types in the order that will be used by the glm", Validate.REQUIRED);
      
    
    @Override
    public void initAndValidate() {
    	String[] splittedTypes = typesInput.get().split("\\s+");

		dimensionInput.set(splittedTypes.length);
    	
		traitToType = new HashMap<>();
		reverseTraitToType = new HashMap<>();
		for (int i = 0; i < splittedTypes.length; i++)
			traitToType.put(splittedTypes[i], i);
		for (int i = 0; i < splittedTypes.length; i++)
			reverseTraitToType.put(i, splittedTypes[i]);
    }

    /**
     * Returns the time to the next interval.
     */
    public double getInterval(int i) {
    	if (i > rateShiftsInput.get().getDimension())
    		return Double.NEGATIVE_INFINITY;
		else
			return rateShiftsInput.get().getArrayValue(i);
    }   

    public boolean intervalIsDirty(int i){
		if (NeGLMInput.get()!=null){
			if(NeGLMInput.get().isDirty())
				return true;
		}else{
			for (int j = 0; j < NeInput.get().getDimension(); j++)
				if (NeInput.get().isDirty(i))
					return true;
		}
		if (migrationGLMInput.get()!=null){
			if(migrationGLMInput.get().isDirty())
				return true;
		}else{
			for (int j = 0; j < migrationRatesInput.get().getDimension(); j++)
				if (migrationRatesInput.get().isDirty(i))
					return true;
		}


    	return false;

    }  
    
    public double arrayMax(double[] arr) {
        double max = Double.NEGATIVE_INFINITY;

        for(double cur: arr)
            max = Math.max(max, cur);

        return max;
    }

    
	@Override
    public double[] getCoalescentRate(int i){
		if (NeGLMInput.get()!=null){
			double[] Ne = NeGLMInput.get().getRates();
			double[] coal = new double[Ne.length];
			for (int j = 0; j < Ne.length; j++)
				coal[j] = 1/Ne[j];
			
			return coal;
		}else{
			double[] coal = new double[dimensionInput.get()];
			for (int j = 0; j < NeInput.get().getDimension(); j++)
				coal[j] = 1/NeInput.get().getArrayValue(j);
				
			return coal;
		}
    }
    
	@Override    
    public double[][] getBackwardsMigration(int i){
    	double[][] m = new double[dimensionInput.get()][dimensionInput.get()];
		if (migrationGLMInput.get()!=null){
			double[] mig = migrationGLMInput.get().getRates();
			int c = 0;
			for (int a = 0; a < dimensionInput.get(); a++){
				for (int b = 0; b < dimensionInput.get(); b++){
					if (a!=b){
						m[b][a] = mig[c];
						c++;
					}
				}
			}
			return m;
		}else{
			int c = 0;
			for (int a = 0; a < dimensionInput.get(); a++){
				for (int b = 0; b < dimensionInput.get(); b++){
					if (a!=b){
						m[b][a] = migrationRatesInput.get().getArrayValue(c);
						c++;
					}
				}
			}
			return m;
		}
    }

	@Override
	public void recalculate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void init(PrintStream out) {
		for (int i = 0 ; i < dimensionInput.get(); i++){
			out.print(String.format("Ne.%s.%s\t", getID(), getStringStateValue(i)));
		}
    	int c = 0;
    	for (int a = 0; a < dimensionInput.get(); a++){
    		for (int b = 0; b < dimensionInput.get(); b++){
    			if (a!=b){
					out.print(String.format("migration.%s.%s_to_%s\t", getID(), getStringStateValue(a), getStringStateValue(b)));
    				c++;
    			}
    		}
    	} 
		
		

	}

	@Override
	public void log(int sample, PrintStream out) {
		double[] Ne = new double[dimensionInput.get()];
		if (NeGLMInput.get()!=null){
			Ne = NeGLMInput.get().getRates();
		}else{
			for (int j = 0; j < NeInput.get().getDimension(); j++)
				Ne[j] = NeInput.get().getArrayValue(j);
		}
		double[][] m = new double[dimensionInput.get()][dimensionInput.get()];
		if (migrationGLMInput.get()!=null){
			double[] mig = migrationGLMInput.get().getRates();
			int c = 0;
			for (int a = 0; a < dimensionInput.get(); a++){
				for (int b = 0; b < dimensionInput.get(); b++){
					if (a!=b){
						m[a][b] = mig[c];
						c++;
					}
				}
			}
		}else{
			int c = 0;
			for (int a = 0; a < dimensionInput.get(); a++){
				for (int b = 0; b < dimensionInput.get(); b++){
					if (a!=b){
						m[a][b] = migrationRatesInput.get().getArrayValue(c);
						c++;
					}
				}
			}
		}



		for (int i = 0 ; i < dimensionInput.get(); i++){
			out.print(Ne[i] + "\t");
		}
		
    	int c = 0;
    	for (int a = 0; a < dimensionInput.get(); a++){
    		for (int b = 0; b < dimensionInput.get(); b++){
    			if (a!=b){
					out.print(m[a][b] + "\t");
    				c++;
    			}
    		}
    	} 
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	

    
}