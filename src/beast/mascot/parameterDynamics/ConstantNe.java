package beast.mascot.parameterDynamics;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.mascot.parameterDynamics.EffectivePopulationSizeDynamics;

public class ConstantNe extends EffectivePopulationSizeDynamics {
	
    public Input<RealParameter> NeInput = new Input<>(
    		"logNe", "input of the Ne at the time of the most recent sampled ancestor", Validate.REQUIRED);    

    RealParameter Ne;
    
	@Override
	public void initAndValidate() {
		// should be called by a time
		isTime = true;
	
		Ne = NeInput.get();
	}

	@Override
	public void recalculate() {
	}

	@Override
	public double getNeTime(double t) {
		return Math.exp(Ne.getArrayValue());
	}
	
	@Override
	public boolean isDirty() {
		if (Ne.isDirty(0))
			return true;
		
		return false;
	}
}
