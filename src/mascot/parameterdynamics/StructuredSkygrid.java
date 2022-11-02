package mascot.parameterdynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.EffectivePopulationSizeDynamics;

@Description("Skygrid style dynamics for MASCOT. These assume constant effective population sizes within one state "+
				"and that these effective population sizes only change within one interval")
public class StructuredSkygrid extends EffectivePopulationSizeDynamics {
	
    public Input<RealParameter> NeLogInput = new Input<>(
    		"NeLog", "input of the log effective population sizes", Validate.REQUIRED);    

    RealParameter NeLog;
    
	@Override
	public void initAndValidate() {
		// should be called by an interval
		isTime = false;
		NeLog = NeLogInput.get();
	}
	
	@Override
	public void setNrIntervals(int intervals) {
		if (NeLogInput.get().getDimension()!=intervals) 
			NeLogInput.get().setDimension(intervals);		
	}

	@Override
	public void recalculate() {
	}
	
	
	public double getNeInterval(int i) {
		// get in which interval the current time falls		
		return Math.exp(NeLog.getArrayValue(i));
	}

	@Override
	public boolean isDirty() {
		for (int i = 0; i < NeLog.getDimension(); i++)
			if(NeLogInput.get().isDirty(i))
				return true;

		return false;	
	}

}
