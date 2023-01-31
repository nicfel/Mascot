package mascot.util;

import beast.base.core.Input;
import mascot.dynamics.RateShifts;


/**
 * Taken over from TIMOTHY VAUGHAN
 */
public class InitializedRateShiftsForSkyline extends RateShifts {
	
    public InitializedRateShiftsForSkyline() {   
    	valuesInput.setRule(Input.Validate.OPTIONAL);
    }

	@Override
	public void initAndValidate() {
		super.initAndValidate();
	}
}
