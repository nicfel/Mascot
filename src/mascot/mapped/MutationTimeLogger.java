package mascot.mapped;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

import java.io.PrintStream;


@Description("Reconstructs sequences at internal nodes and logs them in NEXUS format")
public class MutationTimeLogger extends BEASTObject implements Loggable {

	public Input<AncestralMappedTreeSequenceLogger> ancestralMappedTreeSequenceLoggerInput = new Input<AncestralMappedTreeSequenceLogger>(
			"ancestralMappedTreeSequenceLogger", "AncestralMappedTreeSequenceLogger", Input.Validate.REQUIRED);

	AncestralMappedTreeSequenceLogger amts;


	@Override
	public void initAndValidate() {
		amts = ancestralMappedTreeSequenceLoggerInput.get();
	}

	@Override
	public void init(PrintStream out) {
		out.print("Sample\tNT\tPosition");
		for (int i = 0; i < amts.mappedMascotInput.get().getDimension(); i++){
			out.print("\t" + amts.mappedMascotInput);
		}

	}

	@Override
	public void log(long sample, PrintStream out) {
		amts.reMapTree(sample);

		double[][][] time = new double[amts.getNodeSequence(amts.mappedTree.getRoot().getParent()).length]
				[amts.dataType.getStateCount()][amts.mappedMascotInput.get().getDimension()];

		getBranchMutationLength(amts.mappedTree.getRoot(), time);

		// for each position in the sequence, check if there has been a mutation
		for (int i = 0; i < time.length; i++){
			double[] timeNt = new double[time[i].length];
			for (int j = 0; j < time[i].length; j++){
				for (int k = 0; k < time[i][j].length; k++){
					timeNt[j] += time[i][j][k];
				}
			}

			int mutCount = 0;
			for (int j = 0; j < timeNt.length; j++){
				if (timeNt[j] > 0){
					mutCount++;
				}
			}

			if (mutCount>1){
				for (int j = 0; j < time[i].length; j++){
					out.print(sample);
					if (timeNt[j] > 0){
						for (int k = 0; k < time[i][j].length; k++) {
							out.print(sample + "\t");
						}
					}
				}
			}

		}
	}

	private void getBranchMutationLength(Node node, double[][][] time) {
		if (node.isLeaf()) {
			return;
		}

		if (node.getParent()!=null) {
			int[] parentSequence = amts.getNodeSequence(node.getParent());
			int[] childSequence = amts.getNodeSequence(node);
			for (int i = 0; i < parentSequence.length; i++) {
				if (parentSequence[i] != childSequence[i]) {
					double urand = Randomizer.nextDouble();
					time[i][parentSequence[i]][(int) node.getParent().getMetaData("location")] += node.getLength() * urand;
					time[i][childSequence[i]][(int) node.getParent().getMetaData("location")] += node.getLength() * (1 - urand);
				} else {
					time[i][parentSequence[i]][(int) node.getParent().getMetaData("location")] += node.getLength();
				}
			}
		}

		if (node.getChildCount() == 1) {
			throw new IllegalArgumentException("node " + node.getID() + " has only one child");
		} else {
			getBranchMutationLength(node.getLeft(), time);
			getBranchMutationLength(node.getRight(), time);
		}
	}

	private void appendDouble(StringBuffer buf, double d) {
		buf.append(d);
	}

	@Override
	public void close(PrintStream out) {
	}
}