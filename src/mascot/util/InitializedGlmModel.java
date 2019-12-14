package mascot.util;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TraitSet;
import beast.mascot.glmmodel.GlmModel;
import beast.mascot.glmmodel.LogLinear;

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
public class InitializedGlmModel extends LogLinear {

    public InitializedGlmModel() {    	
    	covariateListInput.setRule(Input.Validate.OPTIONAL);
    	clockInput.setRule(Input.Validate.OPTIONAL);
    	scalerInput.setRule(Input.Validate.OPTIONAL);
    	indicatorInput.setRule(Input.Validate.OPTIONAL);
   }

	@Override
	public void init(PrintStream out) {
		super.init(out);
	}

	@Override
	public void close(PrintStream out) {
		 super.close(out);		
	}

	@Override
	public double[] getRates(int i) {
		return super.getRates(i);
	}

	@Override
	public void initAndValidate() {
		super.initAndValidate();
	}
}
