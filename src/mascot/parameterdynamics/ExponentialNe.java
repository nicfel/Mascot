package mascot.parameterdynamics;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

public class ExponentialNe extends NeDynamics {
	
    public Input<RealParameter> logNeNullInput = new Input<>(
    		"NeNull", "input of the Ne at the time of the most recent sampled ancestor", Validate.REQUIRED);    
    public Input<RealParameter> growthRateInput = new Input<>(
    		"growthRate", "input of the growth rate", Validate.REQUIRED);    

    public Input<Double> minNeInput = new Input<>(
    		"minNe", "input of the minimal Ne", 0.0);    

    RealParameter logNeNull;
    RealParameter growthRate;
    
	@Override
	public void initAndValidate() {
		// should be called by a time
		isTime = true;
		
		logNeNull = logNeNullInput.get();
		growthRate = growthRateInput.get();
	}

	@Override
	public void recalculate() {
	}

	@Override
	public double getNeTime(double t) {
		
		return Math.max(minNeInput.get(), Math.exp(logNeNull.getArrayValue()-t*growthRate.getArrayValue()));
	}
	
	@Override
	public boolean isDirty() {
		if (logNeNull.isDirty(0))
			return true;
		
		if (growthRate.isDirty(0))
			return true;
		
		return false;
	}


}
