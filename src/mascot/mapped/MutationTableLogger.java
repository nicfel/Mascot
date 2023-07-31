package mascot.mapped;

import beast.base.core.*;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.io.PrintStream;


@Description("Reconstructs sequences at internal nodes and logs them in NEXUS format")
public class MutationTableLogger extends BEASTObject implements Loggable {

	public Input<AncestralMappedTreeSequenceLogger> ancestralMappedTreeSequenceLoggerInput = new Input<AncestralMappedTreeSequenceLogger>(
			"ancestralMappedTreeSequenceLogger", "AncestralMappedTreeSequenceLogger", Input.Validate.REQUIRED);


	AncestralMappedTreeSequenceLogger amts;
	boolean first;

	@Override
	public void initAndValidate() {
		amts = ancestralMappedTreeSequenceLoggerInput.get();
	}

	@Override
	public void init(PrintStream out) {
		out.print("Mutation\tAA\tBranchType\tBranchLength\tBranchMutations");
	}

	@Override
	public void log(long sample, PrintStream out) {
		amts.reMapTree(sample);

		Tree mt = amts.mappedTree.copy();
		amts.pruneSingleChildNodes(mt);
		first = true;
		getBranchMutations(sample, out, mt.getRoot());
	}

	private void getBranchMutations(long sample, PrintStream out, Node node) {
		if (node.isLeaf()) {
			return;
		}
		if (node.getParent()!=null) {
			// get the mutations on the branch
			int[] parentSequence = amts.getNodeSequence(node.getParent());
			int[] childSequence = amts.getNodeSequence(node);

			int totMutations = 0;
			for (int i = 0; i < parentSequence.length; i++) {
				if (parentSequence[i] != childSequence[i]) {
					totMutations+=1;
				}
			}
			for (int i = 0; i < parentSequence.length; i++) {
				if (parentSequence[i] != childSequence[i]) {
					if (!first)
						out.println(sample + "\t");
					first = false;
					String mut = amts.getMutation(i, parentSequence[i], childSequence[i]);
					out.print(mut + "\t" +
							amts.getAAchange(i, parentSequence, childSequence) + "\t" + amts.getBranchType(node)  + "\t" + node.getLength() + "\t" + totMutations);
					out.print("\n");
				}
			}
		}


		if (node.getChildCount() == 1) {
			getBranchMutations(sample, out, node.getChild(0));
		} else {
			getBranchMutations(sample, out, node.getLeft());
			getBranchMutations(sample, out, node.getRight());
		}
	}

	private void appendDouble(StringBuffer buf, double d) {
		buf.append(d);
	}

	@Override
	public void close(PrintStream out) {
	}
}