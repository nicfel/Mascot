package beast.mascot.glmmodel;

import java.util.ArrayList;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;

public abstract class GlmModel extends CalculationNode implements Loggable {
	
    public Input<List<RealParameter>> covariatesInput = new Input<>("covariates", "input of covariates", new ArrayList<>(), Validate.REQUIRED);
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
			int l1 = i*(dim*(dim-1));
			int l2 = (i+1)*(dim*(dim-1));
			
			// check that the dimension of the covariates are correct
			for (int j = 0; j < covariatesInput.get().size(); j++){
				if (covariatesInput.get().get(j).getDimension()!=l1 &&
						covariatesInput.get().get(j).getDimension()!=l2)
				throw new RuntimeException("The dimension of the the covariate \"" + covariatesInput.get().get(j).getID() + "\" is wrong.\n" + 
						"The current dimension is " + covariatesInput.get().get(j).getDimension() +
						", but should be equal to the number of rate shifts\n"+
						"(time intervals or time intervals+1) times states(states-1)\n" +
						"i.e. it should either be " + l1 + " or " +l2 + "\n");
			}
		}else{
			// calc the two possible lengths of the covariates
			int l1 = i*dim;
			int l2 = (i+1)*dim;
			for (int j = 0; j < covariatesInput.get().size(); j++){
				if (covariatesInput.get().get(j).getDimension()!=l1 &&
						covariatesInput.get().get(j).getDimension()!=l2)
				throw new RuntimeException("The dimension of the the covariate \"" + covariatesInput.get().get(j).getID() + "\" is wrong.\n" + 
						"The current dimension is " + covariatesInput.get().get(j).getDimension() +
						", but should be equal to the number of rate shifts\n"+
						"(time intervals or time intervals+1) times states\n" +
						"i.e. it should either be " + l1 + " or " +l2 + "\n");
			}

		}
		
		
		verticalEntries = covariatesInput.get().get(0).getDimension()/nrIntervals;
	}
	

}
