package mascot.parameterdynamics;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import mascot.glmmodel.CovariateList;

public class LogLinearGLM extends NeDynamics {	
    public Input<CovariateList> covariateListInput = new Input<>("covariateList", "input of covariates", Validate.REQUIRED);
    public Input<RealParameter> scalerInput = new Input<>("scaler", "input of covariates scaler", Validate.REQUIRED);    
    public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "input of covariates scaler", Validate.REQUIRED);
    public Input<RealParameter> clockInput = new Input<>("clock", "clock rate of the parameter",Validate.REQUIRED);
    public Input<RealParameter> errorInput = new Input<>("error", "time variant error term in the GLM model for the rates");
    public Input<RealParameter> constantErrorInput = new Input<>("constantError", "time invariant error term in the GLM model for the rates");
    
    final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts","When to switch between elements of Ne", Input.Validate.REQUIRED);

    RealParameter Ne;
    RealParameter rateShifts;
    
    boolean valuesKnown = false;
    double[] rates;


	@Override
	public void initAndValidate() {
		// set the dimension of the scalers, indicators and potentially the error term
    	scalerInput.get().setDimension(covariateListInput.get().size());
    	indicatorInput.get().setDimension(covariateListInput.get().size());
    	
    	if (errorInput.get()!=null)
    		errorInput.get().setDimension(covariateListInput.get().get(0).getDimension());
    	
		isTime = true;
    	rateShifts = rateShiftsInput.get();
    	rates = new double[covariateListInput.get().get(0).getDimension()];
	}
	
	
	@Override
	public double getNeTime(double t) {	
//		if (!valuesKnown)
			recalculate();
		
		int intervalnr = getIntervalNr(t);
		if (intervalnr>=rateShifts.getDimension()) {
			return rates[rateShifts.getDimension()-1];
		}				
		return rates[intervalnr];
	}
	
	
	private int getIntervalNr(double t) {
		// check which interval t + offset is in
		for (int i = 0; i < rateShifts.getDimension(); i++)
			if (t<rateShifts.getArrayValue(i))
				return i;
		
		// after the last interval, just keep using the last element
		return rateShifts.getDimension();					
	}



	@Override
	public void recalculate() {
		for (int i = 0; i < rates.length; i++) {
			double logrates = 0;

			for (int j = 0; j < covariateListInput.get().size(); j++){
				if (indicatorInput.get().getArrayValue(j) > 0.0){
					logrates += scalerInput.get().getArrayValue(j)
							*covariateListInput.get().get(j).getArrayValue(i);				
				}
			}
		
	    	if (errorInput.get()!=null)
	   			logrates += errorInput.get().getArrayValue(i);    	
	    	
	    	if (constantErrorInput.get()!=null)
	   			logrates += constantErrorInput.get().getArrayValue();
	    	
	    	rates[i] = clockInput.get().getArrayValue()*Math.exp(logrates);
		}		
		
		valuesKnown = true;
	}


	@Override
	public boolean isDirty() {
		for (int i = 0; i < scalerInput.get().getDimension(); i++)
			if(scalerInput.get().isDirty(i)){
				valuesKnown = false;
				return true;
			}
		
		for (int i = 0; i < indicatorInput.get().getDimension(); i++)
			if(indicatorInput.get().isDirty(i)){
				valuesKnown = false;
				return true;
			}
		
		if (errorInput.get() != null)
			for (int i = 0; i < errorInput.get().getDimension(); i++)
				if(errorInput.get().isDirty(i)){
					valuesKnown = false;
					return true;
				}
		
		if (constantErrorInput.get() != null)
			for (int i = 0; i < constantErrorInput.get().getDimension(); i++)
				if(constantErrorInput.get().isDirty(i)) {
					valuesKnown = false;
					return true;
				}

		
		if (clockInput.get().isDirty(0)) {
			valuesKnown = false;
			return true;
		}
		
		return false;
	}

	@Override
	public void restore() {
		valuesKnown = false;
	}


}
