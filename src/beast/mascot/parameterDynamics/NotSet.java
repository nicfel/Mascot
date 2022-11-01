package beast.mascot.parameterDynamics;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.mascot.parameterDynamics.EffectivePopulationSizeDynamics;

public class NotSet extends EffectivePopulationSizeDynamics {
	
    public Input<Boolean> notSet = new Input<>(
    		"button that does nothing", "notSet", false);    

	
	@Override
	public void initAndValidate() {
		// should be called by a time
		isTime = true;
	}
	
	@Override
	public void recalculate() {
	}

	@Override
	public double getNeTime(double t) {
		return Math.exp(-1);
	}
	
	@Override
	public boolean isDirty() {
		return false;
	}

}
