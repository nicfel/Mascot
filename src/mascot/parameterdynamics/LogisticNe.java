package mascot.parameterdynamics;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

public class LogisticNe extends NeDynamics {
	
    public Input<RealParameter> carryingProportionInput = new Input<>(
    		"carryingProportion", "the proportion of the current Ne of the maximial Ne (capactity)", Validate.REQUIRED);    
    public Input<RealParameter> capacityInput = new Input<>(
    		"capacity", "input of the maximal Ne", Validate.REQUIRED);    
    public Input<RealParameter> growthRateInput = new Input<>(
    		"growthRate", "input of the growth rate", Validate.REQUIRED);  
    
    public Input<Double> minNeInput = new Input<>(
    		"minNe", "input of the minimal Ne", 0.0);    

    RealParameter cP;
    RealParameter capacity;
    RealParameter growthRate;
    
	@Override
	public void initAndValidate() {
		// should be called by a time
		isTime = true;
		
		cP = carryingProportionInput.get();
		
		capacity = capacityInput.get();
		growthRate = growthRateInput.get();		
	}

	@Override
	public void recalculate() {
	}

	@Override
	public double getNeTime(double t) {		
		return Math.max(minNeInput.get(),  Math.exp(capacity.getValue())/(1 + (1-cP.getValue())/cP.getValue() * Math.exp(t*growthRate.getArrayValue())));
	}

	@Override
	public boolean isDirty() {
		if (cP.isDirty(0))
			return true;
		
		if (capacity.isDirty(0))
			return true;

		if (growthRate.isDirty(0))
			return true;
		
		return false;
	}	

}
