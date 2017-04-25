
package beast.mascot.distribution;


import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.evolution.tree.TreeInterface;

// should be the same as the original tree distribution but allowing for Structured tree intervals 
// additionally allows recalculation of which daughter lineages were involved in a coalescent event
@Description("Distribution on a tree, typically a prior such as Coalescent or Yule")
public class StructuredTreeDistribution extends Distribution {
    public Input<TreeInterface> treeInput = new Input<TreeInterface>("tree", "tree over which to calculate a prior or likelihood");
    public Input<StructuredTreeIntervals> treeIntervalsInput = new Input<StructuredTreeIntervals>("structuredTreeIntervals",
    		"Structured Intervals for a phylogenetic beast tree", Validate.XOR, treeInput);

    
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
        final StructuredTreeIntervals ti = treeIntervalsInput.get();
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
