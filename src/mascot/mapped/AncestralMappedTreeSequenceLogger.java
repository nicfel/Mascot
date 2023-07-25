package mascot.mapped;

import beast.base.core.*;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.FilteredAlignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import mascot.distribution.MappedMascot;

import java.io.PrintStream;
import java.util.*;


@Description("Reconstructs sequences at internal nodes and logs them in NEXUS format")
public class AncestralMappedTreeSequenceLogger extends BEASTObject implements Function, Loggable {

	public Input<String> tagInput = new Input<String>("tag", "label used to report trait", "mutations");

	final public Input<Boolean> logMutationsOnlyInput = new Input<>("logMutationsOnly",
			"if true, logs only mutations on edges instead of sequences on nodes", false);

	final public Input<List<AncestralStateTreeLikelihood>> ancestralTreeLikelihoodInput = new Input<>(
			"ancestralTreeLikelihood",
			"ancestral tree likelihoods that," + "if more than one is inputted is assumed to be a filtered alignment",
			new ArrayList<>());

	final public Input<Boolean> translateInput = new Input<>("translate", "if true, dna is translated to amino acids",
			false);

	public Input<MappedMascot> mappedMascotInput = new Input<MappedMascot>("mappedMascot", "list of leaf traits",
			Input.Validate.REQUIRED);

//	int [] siteStates;

	protected DataType dataType;

	int totalLength;
	List<Integer[]> mapping;

	Map<String, String> translationsTable;

	@Override
	public void init(PrintStream out) {
		((Tree) ancestralTreeLikelihoodInput.get().get(0).treeInput.get()).init(out);

		dataType = ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getDataType();
		totalLength = ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getSiteCount();
		;

		mapping = new ArrayList<>();

		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			if (ancestralTreeLikelihoodInput.get().get(i).dataInput.get() instanceof FilteredAlignment) {
				mapping.add(parseFilter(
						((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).filterInput
								.get(),
						((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).alignmentInput
								.get(),
						ancestralTreeLikelihoodInput.get().get(i).dataInput.get().getSiteCount()));

				totalLength = ((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput
						.get()).alignmentInput.get().getSiteCount();
			} else {
				Integer[] newMap = new Integer[totalLength];
				for (int j = 0; j < newMap.length; j++)
					newMap[j] = j;
				mapping.add(newMap);
			}
		}
		if (translateInput.get())
			initTable();
	}

	private Integer[] parseFilter(String filterString, Alignment alignment, int redCount) {
		// parse filter specification
//        String filterString = filterInput.get();
		String[] filters = filterString.split(",");
		int[] from = new int[filters.length];
		int[] to = new int[filters.length];
		int[] step = new int[filters.length];
		for (int i = 0; i < filters.length; i++) {
			filterString = " " + filters[i] + " ";
			if (filterString.matches(".*-.*")) {
				// range, e.g. 1-100/3
				if (filterString.indexOf('\\') >= 0) {
					String str2 = filterString.substring(filterString.indexOf('\\') + 1);
					step[i] = parseInt(str2, 1);
					filterString = filterString.substring(0, filterString.indexOf('\\'));
				} else {
					step[i] = 1;
				}
				String[] strs = filterString.split("-");
				from[i] = parseInt(strs[0], 1) - 1;
				to[i] = parseInt(strs[1], alignment.getSiteCount()) - 1;
			} else if (filterString.matches(".*:.*:.+")) {
				// iterator, e.g. 1:100:3
				String[] strs = filterString.split(":");
				from[i] = parseInt(strs[0], 1) - 1;
				to[i] = parseInt(strs[1], alignment.getSiteCount()) - 1;
				step[i] = parseInt(strs[2], 1);
			} else if (filterString.trim().matches("[0-9]*")) {
				from[i] = parseInt(filterString.trim(), 1) - 1;
				to[i] = from[i];
				step[i] = 1;
			} else {
				throw new IllegalArgumentException("Don't know how to parse filter " + filterString);
			}
		}

		boolean[] used = new boolean[alignment.getSiteCount()];
		for (int i = 0; i < to.length; i++) {
			for (int k = from[i]; k <= to[i]; k += step[i]) {
				used[k] = true;
			}
		}
		Integer[] map = new Integer[redCount];
		int c = 0;
		for (int i = 0; i < used.length; i++) {
			if (used[i]) {
				map[c] = i;
				c++;
			}
		}

		return map;
	}

	int parseInt(String str, int defaultValue) {
		str = str.replaceAll("\\s+", "");
		try {
			return Integer.parseInt(str);
		} catch (Exception e) {
			return defaultValue;
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			ancestralTreeLikelihoodInput.get().get(i).calculateLogP();
			ancestralTreeLikelihoodInput.get().get(i).redrawAncestralStates();
		}
		mappedMascotInput.get().calculateLogP();
		Node mappedRoot = mappedMascotInput.get().getRoot();
		mappedRoot.sort();
		// compute mutations for mapped tree
		mapMuts(mappedRoot, null);

		out.print("tree STATE_" + sample + " = ");

		out.print(toNewick(mappedRoot, null, null));
		out.print(";");
	}

	private void mapMuts(Node node, int[] parentSiteStates) {
		int[] currSiteStates = new int[totalLength];

		for (int i = 0; i < currSiteStates.length; i++)
			currSiteStates[i] = -1;

		// If the node has only one child, it is a migration event
		if (node.getChildCount() == 1) {
			List<Double> times = new ArrayList<>();

			// gather child times until the next coal or leaf node, go to the next
			double currtime = node.getParent().getHeight();
			Node currNode = node;
			while (currNode.getChildCount() == 1) {
				times.add(currtime - currNode.getHeight());
				currtime = currNode.getHeight();
				currNode = currNode.getChild(0);
			}
			times.add(currtime - currNode.getHeight());

			// get the sequence at the currtime
			for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
				int[] siteStates = ancestralTreeLikelihoodInput.get().get(i)
						.getStatesForNode(ancestralTreeLikelihoodInput.get().get(0).treeInput.get(), currNode);

				for (int j = 0; j < mapping.get(i).length; j++) {
					currSiteStates[mapping.get(i)[j]] = siteStates[j];
				}
			}

			// get all the mutations and assign them to intervals
			double[] array = new double[times.size()];
			for (int i = 0; i < times.size(); i++)
				array[i] = times.get(i);

			// look for differences
			List<String> muts = new ArrayList<>();
			List<Integer> ints = new ArrayList<>();
			for (int i = 0; i < parentSiteStates.length; i++) {
				if (parentSiteStates[i] != currSiteStates[i]) {
					muts.add(dataType.getCharacter(parentSiteStates[i]) + "" + (i + 1) + ""
							+ dataType.getCharacter(currSiteStates[i]));
					ints.add(Randomizer.randomChoicePDF(array));
				}
			}

			int k = 0;
			int l = 0;
			while (node.getChildCount() == 1) {
				String muts_string = "\"";
				boolean print=false;
				for (int i = 0; i < ints.size(); i++) {
					if (ints.get(i) == k) {
						muts_string = muts_string + "," + muts.get(i);
						print = true;
						l++;
					}
				}
				muts_string = muts_string + "\"";
				muts_string = muts_string.replace("\",", "\"");
				if (print)
					node.setMetaData(tagInput.get(), muts_string);
				k++;
				node = node.getChild(0);
			}
			
			String muts_string = "\"";
			boolean print=false;

			for (int i = 0; i < ints.size(); i++) {
				if (ints.get(i) == k) {
					muts_string = muts_string + "," + muts.get(i);
					print = true;
					l++;
				}
			}
			muts_string = muts_string + "\"";
			muts_string = muts_string.replace("\",", "\"");
			if (print)
				node.setMetaData(tagInput.get(), muts_string);

		} else {
			// node is leaf or coalescent event
			for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
				int[] siteStates = ancestralTreeLikelihoodInput.get().get(i)
						.getStatesForNode(ancestralTreeLikelihoodInput.get().get(0).treeInput.get(), node);

				for (int j = 0; j < mapping.get(i).length; j++) {
					currSiteStates[mapping.get(i)[j]] = siteStates[j];
				}
			}

			if (!node.isRoot()) {
				List<String> muts = new ArrayList<>();
				for (int i = 0; i < parentSiteStates.length; i++) {
					if (parentSiteStates[i] != currSiteStates[i]) {
						muts.add(dataType.getCharacter(parentSiteStates[i]) + "" + (i + 1) + ""
								+ dataType.getCharacter(currSiteStates[i]));
					}
				}

				String muts_string = "\"";
				boolean print=false;

				for (int i = 0; i < muts.size(); i++) {
					print = true;
					muts_string = muts_string + "," + muts.get(i);
				}
				muts_string = muts_string + "\"";
				muts_string=muts_string.replace("\",", "\"");
				if (print)
					node.setMetaData(tagInput.get(), muts_string);
			}
		}

		if (!node.isLeaf()) {
			mapMuts(node.getLeft(), currSiteStates);
			mapMuts(node.getRight(), currSiteStates);
		}
	}

	private void appendDouble(StringBuffer buf, double d) {
		buf.append(d);
	}

	String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
		StringBuffer buf = new StringBuffer();
		if (node.getLeft() != null) {
			buf.append("(");
			buf.append(toNewick(node.getLeft(), metadataList, branchRateModel));
			if (node.getRight() != null) {
				buf.append(',');
				buf.append(toNewick(node.getRight(), metadataList, branchRateModel));
			}
			buf.append(")");
		} else {
			buf.append(node.getNr() + 1);
		}
		if (!node.isLeaf()) {

			buf.append("[&");
			buf.append(mappedMascotInput.get().dynamics.typeTraitInput.getName() + "="
					+ mappedMascotInput.get().dynamics.getStringStateValue((int) node.getMetaData("location")));

			if (node.getMetaData(tagInput.get())!=null)
				buf.append("," + tagInput.get() + "=" + node.getMetaData(tagInput.get()));
			

			if (branchRateModel != null) {
				buf.append(",rate=");
				appendDouble(buf, branchRateModel.getRateForBranch(node));
			}
			buf.append(']');

		} else {
			String sampleID = node.getID();
			String[] splits = sampleID.split("_");
			int sampleState;

			if (mappedMascotInput.get().dynamics.typeTraitInput.get() != null) {
				sampleState = mappedMascotInput.get().dynamics.getValue(node.getID());
			}

			else {
				sampleState = Integer.parseInt(splits[splits.length - 1]); // samples states (or priors) should
																			// eventually be specified in the XML
			}

//			if (node.getMetaData("location") != null) {
				buf.append("[&");
				buf.append(mappedMascotInput.get().dynamics.typeTraitInput.getName() + "="
						+ mappedMascotInput.get().dynamics.getStringStateValue((int) node.getMetaData("location")));
				if (node.getMetaData(tagInput.get())!=null)
					buf.append("," + tagInput.get() + "=" + node.getMetaData(tagInput.get()));
				
				
				buf.append(']');
//			}

		}

		buf.append(":");
//		if (substitutions) {
//			appendDouble(buf, node.getLength() * branchRateModel.getRateForBranch(node));
//		} else {
		appendDouble(buf, node.getLength());
//		}
		return buf.toString();
	}

	@Override
	public void close(PrintStream out) {
		((Tree) ancestralTreeLikelihoodInput.get().get(0).treeInput.get()).close(out);
	}

	@Override
	public int getDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getArrayValue(int dim) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub

	}

	private void initTable() {
		String dna = "AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT";
		String aa = "K N K N T T T T R S R S I I M I Q H Q H P P P P R R R R L L L L E D E D A A A A G G G G V V V V O Y O Y S S S S O C W C L F L F";

		translationsTable = new HashMap<>();
		String[] dnaarray = dna.split(" ");
		String[] aaarray = aa.split(" ");
		for (int i = 0; i < dnaarray.length; i++)
			translationsTable.put(dnaarray[i], aaarray[i]);

	}
}