
package mascot.distribution;


import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeIntervals;

// should be the same as the original tree distribution but allowing for Structured tree intervals 
// additionally allows recalculation of which daughter lineages were involved in a coalescent event
@Description("Distribution on a tree, typically a prior such as Coalescent or Yule")
public class StructuredTreeDistribution extends TreeDistribution {

    public Input<StructuredTreeIntervals> structuredTreeIntervalsInput = new Input<StructuredTreeIntervals>("structuredTreeIntervals",
    		"Structured Intervals for a phylogenetic beast tree", Validate.REQUIRED);
    
 
    
    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    protected boolean requiresRecalculation() {
        final StructuredTreeIntervals ti = structuredTreeIntervalsInput.get();
        if (ti != null) {
            assert ti.isDirtyCalculation();
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }
    
 	/** Indicate that the tree distribution can deal with dated tips in the tree
	 * Some tree distributions like the Yule prior cannot handle this.
	 * @return true by default
	 */
	public boolean canHandleTipDates() {
		return true;
	}
}
