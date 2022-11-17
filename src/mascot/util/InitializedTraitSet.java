package mascot.util;

import beast.base.core.Input;
import beast.base.evolution.tree.TraitSet;

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
public class InitializedTraitSet extends TraitSet {

    public InitializedTraitSet() {
        traitsInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        if (traitsInput.get() == null) {
            String value = taxaInput.get().getTaxaNames().stream()
                    .map(n -> n + "=NOT_SET")
                    .collect(Collectors.joining(","));

            traitsInput.setValue(value, this);
        }

        super.initAndValidate();
    }
}
