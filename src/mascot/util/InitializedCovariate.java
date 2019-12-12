package mascot.util;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TraitSet;
import beast.mascot.glmmodel.Covariate;
import beast.mascot.glmmodel.CovariateList;
import beast.mascot.glmmodel.GlmModel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

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
public class InitializedCovariate extends Covariate {

    public InitializedCovariate() {    	
    	valuesInput.setRule(Input.Validate.OPTIONAL);
   }

	@Override
	public void initAndValidate() {
	}
}
