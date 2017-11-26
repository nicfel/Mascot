package beast.mascot.glmmodel;

import java.util.ArrayList;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;

public abstract class GLMStepwiseModel extends CalculationNode implements Loggable {
	
    public Input<List<RealParameter>> covariatesInput = new Input<>("covariates", "input of covariates", new ArrayList<>(), Validate.REQUIRED);
    public Input<RealParameter> scalerInput = new Input<>("scaler", "input of covariates scaler", Validate.REQUIRED);    
    public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "input of covariates scaler", Validate.REQUIRED);
    public Input<RealParameter> clockInput = new Input<>("clock", "clock rate of the parameter",Validate.REQUIRED);
    public Input<RealParameter> errorInput = new Input<>("error", "error term in the GLM model for the- rates");
    
    public int nrIntervals;
    public int verticalEntries;
    
	public abstract double[] getRates(int i);
	
	public boolean isDirty(){
		for (int i = 0; i < scalerInput.get().getDimension(); i++)
			if(scalerInput.get().isDirty(i) && indicatorInput.get().isDirty(i))
					return true;
		
		for (int i = 0; i < indicatorInput.get().getDimension(); i++)
			if(indicatorInput.get().isDirty(i))
					return true;
				
		if (clockInput.get().isDirty(0))
			return true;
		
		return false;
	}

	public void setNrIntervals(int i){
		nrIntervals = i;
		verticalEntries = covariatesInput.get().get(0).getDimension()/nrIntervals;
	}
	

}
