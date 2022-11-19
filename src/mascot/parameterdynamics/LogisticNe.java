package mascot.parameterdynamics;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

public class LogisticNe extends NeDynamics {
	
    public Input<RealParameter> carryingPercentageInput = new Input<>(
    		"carryingPercentage", "input of the Ne at the time of the most recent sampled ancestor");    
    public Input<RealParameter> logitCarryingPercentageInput = new Input<>(
    		"logitCarryingPercentage", "input of the Ne at the time of the most recent sampled ancestor", Validate.XOR, carryingPercentageInput);    
    public Input<RealParameter> capacityInput = new Input<>(
    		"capacity", "input of the Ne at the time of the most recent sampled ancestor", Validate.REQUIRED);    
    public Input<RealParameter> growthRateInput = new Input<>(
    		"growthRate", "input of the growth rate", Validate.REQUIRED);    


    RealParameter cP;
    RealParameter lCP;
    RealParameter capacity;
    RealParameter growthRate;
    boolean logitInput;
    
	@Override
	public void initAndValidate() {
		// should be called by a time
		isTime = true;
		// check if the final Ne as a percentage of the carrying capacity is logit transformed or not
		if (carryingPercentageInput.get()!=null) {
			cP = carryingPercentageInput.get();
			logitInput = false;
		}else {
			lCP = logitCarryingPercentageInput.get();
			logitInput = true;
		}
		
		capacity = capacityInput.get();
		growthRate = growthRateInput.get();
		
	}

	@Override
	public void recalculate() {
	}

	@Override
	public double getNeTime(double t) {		
		if (logitInput)
			return Math.exp(capacity.getValue())/(1 + (1-getLogitReversed(lCP.getValue()))/getLogitReversed(lCP.getValue()) * Math.exp(t*growthRate.getArrayValue()));
		else
			return Math.exp(capacity.getValue())/(1 + (1-cP.getValue())/cP.getValue() * Math.exp(t*growthRate.getArrayValue()));
	}

	@Override
	public boolean isDirty() {
		if (logitInput)
			if (lCP.isDirty(0))
				return true;
		else
			if (cP.isDirty(0))
				return true;
		
		if (capacity.isDirty(0))
			return true;

		if (growthRate.isDirty(0))
			return true;
		
		return false;
	}
	
	private double getLogitReversed(double logitval) {
		double eLV = Math.exp(logitval);
		return eLV/(1+eLV);
	}

}
