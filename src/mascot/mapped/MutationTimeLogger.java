package mascot.mapped;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;


@Description("Reconstructs sequences at internal nodes and logs them in NEXUS format")
public class MutationTimeLogger extends BEASTObject implements Loggable {

	public Input<AncestralMappedTreeSequenceLogger> ancestralMappedTreeSequenceLoggerInput = new Input<AncestralMappedTreeSequenceLogger>(
			"ancestralMappedTreeSequenceLogger", "AncestralMappedTreeSequenceLogger", Input.Validate.REQUIRED);

	public Input<Boolean> translateInput = new Input<Boolean>("translate", "translate to amino acids", false);

	AncestralMappedTreeSequenceLogger amts;


	@Override
	public void initAndValidate() {
		amts = ancestralMappedTreeSequenceLoggerInput.get();
	}

	@Override
	public void init(PrintStream out) {
		out.print("NT\tPosition");
		for (int i = 0; i < amts.mappedMascotInput.get().dynamicsInput.get().getDimension(); i++){
			out.print("\t" + amts.mappedMascotInput.get().dynamicsInput.get().getStringStateValue(i));
		}

	}

	@Override
	public void log(long sample, PrintStream out) {


		boolean first = true;
		amts.reMapTree(sample);


		double[][][] time = new double[amts.getNodeSequence(amts.mappedTree.getRoot()).length]
				[amts.dataType.getStateCount()][amts.mappedMascotInput.get().dynamicsInput.get().getDimension()];
		if (translateInput.get())
			time = new double[amts.getNodeSequence(amts.mappedTree.getRoot()).length]
					[21][amts.mappedMascotInput.get().dynamicsInput.get().getDimension()];

		getBranchMutationLength(amts.mappedTree.getRoot(), time);

		// for each position in the sequence, check if there has been a mutation
		for (int i = 0; i < time.length; i++){
			// compute the time position i has spent in each nucleotide
			double[] timeNt = new double[time[i].length];
			for (int j = 0; j < time[i].length; j++){
				for (int k = 0; k < time[i][j].length; k++){
					timeNt[j] += time[i][j][k];
				}
			}

			// count the number of mutations based on if position i has spent time in multiple nucleotides
			int mutCount = 0;
			for (int j = 0; j < timeNt.length; j++){
				if (timeNt[j] > 0){
					mutCount++;
				}
			}

			// if there was a mutation, output the site to file
			if (mutCount>1){
				for (int j = 0; j < time[i].length; j++){
					if (timeNt[j] > 0){
						if (!first)
							out.print(sample + "\t");

						first = false;
						if (translateInput.get())
							out.print(amts.aaList.get(j) + "\t" + (i+1) + "\t");
						else
							out.print(amts.dataType.getCharacter(j) + "\t" + (i+1) + "\t");

						for (int k = 0; k < time[i][j].length; k++) {
							out.print(time[i][j][k] + "\t");
						}
						out.print("\n");
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
			// get sequence of parent node
			int[] parentSequence = amts.getNodeSequence(node.getParent());
			if (translateInput.get()){
				parentSequence = amts.getTranslatedNodeSequence(node.getParent());
			}

			List<Node> nodeList = new ArrayList<>();
			// add node to list
			nodeList.add(node);
			// add all single child nodes to list and keep track of total branch length
			double totalBranchLength = getSingleChildList(nodeList);
			// get sequence of child
			int[] childSequence = amts.getNodeSequence(nodeList.get(nodeList.size()-1));
			if (translateInput.get()){
				childSequence = amts.getTranslatedNodeSequence(nodeList.get(nodeList.size()-1));
			}


//			for (Node n : nodeList){
//				System.out.print(n.getChildCount() + " " + n.getLength() + " ");
//			}
//			System.out.print(totalBranchLength + "\n");

			for (int i = 0; i < parentSequence.length; i++) {
				if (parentSequence[i] != childSequence[i]) {
					double mutationTime = Randomizer.nextDouble()*totalBranchLength;
					double currentTime = 0.0;
					for (Node n : nodeList){
						double timeToNextNode = n.getLength();
						if (mutationTime < timeToNextNode){
							time[i][parentSequence[i]][(int) n.getParent().getMetaData("location")] += mutationTime;
							time[i][childSequence[i]][(int) n.getParent().getMetaData("location")] += timeToNextNode-mutationTime;
							timeToNextNode-=mutationTime;
							mutationTime = Double.POSITIVE_INFINITY;
						}else{
							time[i][parentSequence[i]][(int) n.getParent().getMetaData("location")] += timeToNextNode;
							currentTime += timeToNextNode;
							mutationTime -= timeToNextNode;
						}
					}
				} else {
					for (Node n : nodeList)
						time[i][parentSequence[i]][(int) node.getParent().getMetaData("location")] += n.getLength();
				}
			}
			node = nodeList.get(nodeList.size()-1);
		}

		if (node.getChildCount() == 1) {
			throw new IllegalArgumentException("node " + node.getID() + " has only one child");
		} else if (node.getChildCount() == 2) {
			getBranchMutationLength(node.getLeft(), time);
			getBranchMutationLength(node.getRight(), time);
		}
	}

	private double getSingleChildList(List<Node> nodeList) {
		double time = nodeList.get(nodeList.size()-1).getLength();
		if (nodeList.get(nodeList.size()-1).getChildCount()==1) {
			nodeList.add(nodeList.get(nodeList.size() - 1).getChild(0));
			time += getSingleChildList(nodeList);
		}
		return time;
	}

	private void appendDouble(StringBuffer buf, double d) {
		buf.append(d);
	}

	@Override
	public void close(PrintStream out) {
	}
}