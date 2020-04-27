package beast.mascot.glmmodel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;

public abstract class GlmModel extends CalculationNode implements Loggable {
	
    public Input<CovariateList> covariateListInput = new Input<>("covariateList", "input of covariates", Validate.REQUIRED);
    public Input<RealParameter> scalerInput = new Input<>("scaler", "input of covariates scaler", Validate.REQUIRED);    
    public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "input of covariates scaler", Validate.REQUIRED);
    public Input<RealParameter> clockInput = new Input<>("clock", "clock rate of the parameter",Validate.REQUIRED);
    public Input<RealParameter> errorInput = new Input<>("error", "time variant error term in the GLM model for the rates");
    public Input<RealParameter> constantErrorInput = new Input<>("constantError", "time invariant error term in the GLM model for the rates");
    
    public int nrIntervals;
    public int verticalEntries;

    
	public abstract double[] getRates(int i);
	
	public boolean isDirty(){
		for (int i = 0; i < scalerInput.get().getDimension(); i++)
			if(scalerInput.get().isDirty(i))
					return true;
		
		for (int i = 0; i < indicatorInput.get().getDimension(); i++)
			if(indicatorInput.get().isDirty(i))
					return true;
		
		if (errorInput.get() != null)
			for (int i = 0; i < errorInput.get().getDimension(); i++)
				if(errorInput.get().isDirty(i))
						return true;
		
		if (constantErrorInput.get() != null)
			for (int i = 0; i < constantErrorInput.get().getDimension(); i++)
				if(constantErrorInput.get().isDirty(i))
						return true;

		
		if (clockInput.get().isDirty(0))
			return true;
		
		return false;
	}

	public void setNrIntervals(int i, int dim, boolean isMigration){
		nrIntervals = i;
		if (isMigration){
			// calc the two possible lengths of the covariates
			int l2 = i*(dim*(dim-1));
			
			// check that the dimension of the covariates are correct
			for (int j = 0; j < covariateListInput.get().size(); j++){
				if (covariateListInput.get().get(j).getDimension()!=l2)
				throw new RuntimeException("The dimension of the the covariate \"" + covariateListInput.get().get(j).getID() + "\" is wrong.\n" + 
						"The current dimension is " + covariateListInput.get().get(j).getDimension() +
						", but should be equal to the number of rate shifts\n"+
						"i.e. it should be " +l2 + "\n");
			}
		}else{
			// calc the two possible lengths of the covariates
			int l2 = i*dim;
			for (int j = 0; j < covariateListInput.get().size(); j++){
				if (covariateListInput.get().get(j).getDimension()!=l2)
				throw new RuntimeException("The dimension of the the covariate \"" + covariateListInput.get().get(j).getID() + "\" is wrong.\n" + 
						"The current dimension is " + covariateListInput.get().get(j).getDimension() +
						", but should be equal to the number of rate shifts\n"+
						"i.e. it should be " +l2 + "\n");
			}
		}
		
		verticalEntries = covariateListInput.get().get(0).getDimension()/(nrIntervals);
		
	}
	
	public void setNrDummy(){	
		verticalEntries = 0;
	}

	

}
