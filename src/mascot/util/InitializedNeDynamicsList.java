package mascot.util;

import beast.base.core.Input;
import mascot.parameterdynamics.NeDynamicsList;

/**
 * Taken over from TIMOTHY VAUGHAN
 * TraitSet in which value input is optional and the values
 * are initialized to a stand-in value.  Used by the BDMM-Prime BEAUti template,
 * where the trait set must be specified before any value (or indeed
 * the taxa themselves) can be known.
 *
 * Use of this trait set class in the BEAUti template is also important as
 * it causes BDMM-Prime's own traitset input editor to be used rather than
 * some other input editor (such as the one that MTT provdes).
 */
public class InitializedNeDynamicsList extends NeDynamicsList {
	
    public InitializedNeDynamicsList() {    	
        neDynamicsInput.setRule(Input.Validate.OPTIONAL);
   }

	@Override
	public void initAndValidate() {
		super.initAndValidate();
	}
}
