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
import java.util.stream.Collectors;


@Description("Reconstructs sequences at internal nodes and logs them in NEXUS format")
public class AncestralMappedTreeSequenceLogger extends BEASTObject implements Function, Loggable {

	public Input<String> tagInput = new Input<String>("tag", "label used to report trait", "mutations");

	final public Input<Boolean> logMutationsOnlyInput = new Input<>("logMutationsOnly",
			"if true, logs only mutations on edges instead of sequences on nodes", false);

	final public Input<List<AncestralStateTreeLikelihood>> ancestralTreeLikelihoodInput = new Input<>(
			"ancestralTreeLikelihood",
			"ancestral tree likelihoods that," + "if more than one is inputted is assumed to be a filtered alignment",
			new ArrayList<>());

	final public Input<Integer> startReadingFrame = new Input<>("startReadingFrame", "start position of the reading frame for translation");

	final public Input<Integer> endReadingFrame = new Input<>("endReadingFrame", "end position of the reading frame" +
			" if start is specified but not end, stop codons on the root are used for inference", Integer.MAX_VALUE);

	public Input<MappedMascot> mappedMascotInput = new Input<MappedMascot>("mappedMascot", "list of leaf traits",
			Input.Validate.REQUIRED);

	public Input<Boolean> internalNodesOnlyInput = new Input<Boolean>("internalNodesOnly",
			"If true, migration events are not logged default false",
			true);


//	int [] siteStates;

	protected DataType dataType;

	int totalLength;
	List<Integer[]> mapping;

	Map<String, String> translationsTable;
	String stopCodons;
	Tree mappedTree;

	long lastSample = -1;

	String dna = "AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT";
	String aa = "K N K N T T T T R S R S I I M I Q H Q H P P P P R R R R L L L L E D E D A A A A G G G G V V V V O Y O Y S S S S O C W C L F L F";
	List<String> aaList;

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
		if (startReadingFrame.get()!=null)
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

		reMapTree(sample);
		Tree mt = mappedTree.copy();
		pruneSingleChildNodes(mt);

		Node mappedRoot = mt.getRoot();

		mappedRoot.sort();
		// compute mutations for mapped tree
		mapMuts(mappedRoot, null);

		out.print("tree STATE_" + sample + " = ");

		out.print(toNewick(mappedRoot, null, null));
		out.print(";");
	}

	public void reMapTree(long sample){
		if (sample != lastSample) {
			lastSample = sample;
			mappedTree = mappedMascotInput.get().mappedTree.copy();

			for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
				ancestralTreeLikelihoodInput.get().get(i).calculateLogP();
				ancestralTreeLikelihoodInput.get().get(i).redrawAncestralStates();
			}
			mappedMascotInput.get().calculateLogP();
		}
	}

	public void pruneSingleChildNodes(Tree tree){
		if (internalNodesOnlyInput.get()){
			pruneSingleChildNodes(tree.getRoot());
			singleChildCheck(tree.getRoot());
		}
	}

	//  method to check that there are no singlechild nodes left in the tree
	private void singleChildCheck(Node n){
		if (n.getChildCount()==2){
			singleChildCheck(n.getChild(0));
			singleChildCheck(n.getChild(1));
		}else if (n.getChildCount()==1) {
			throw new IllegalArgumentException("single child node left");
		}
	}

	private void pruneSingleChildNodes(Node node) {
		// if the node has only one child, pass down to the next and next until it has more than one child or 0 children
		if (node.getChildCount()==2){
			// keeps track of which child is left or right as that may be changed by the pruning
			Node leftChild = node.getChild(0);
			Node rightChild = node.getChild(1);
			pruneSingleChildNodes(leftChild);
			pruneSingleChildNodes(rightChild);
		} else if (node.getChildCount()==1){
			Node nextNode = getNextNonSingleChildNode(node);
			Node parent = node.getParent();
			nextNode.setParent(parent);
			parent.removeChild(node);
			parent.addChild(nextNode);
			pruneSingleChildNodes(nextNode);
		}
	}

	protected Node getNextNonSingleChildNode(Node node){
		if (node.getChildCount()==1) {
			return getNextNonSingleChildNode(node.getChild(0));
		} else {
			return node;
		}
	}

	private void mapMuts(Node node, int[] parentSiteStates) {
		if (node==null)
			return;

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
//			for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
//				int[] siteStates = ancestralTreeLikelihoodInput.get().get(i)
//						.getStatesForNode(ancestralTreeLikelihoodInput.get().get(0).treeInput.get(), currNode);
//
//				for (int j = 0; j < mapping.get(i).length; j++) {
//					currSiteStates[mapping.get(i)[j]] = siteStates[j];
//				}
//			}
			currSiteStates = getNodeSequence(node);

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

			if (internalNodesOnlyInput.get())
				System.err.println("internalNodesOnlyInput not null, but single child nodes are left at height "+node.getHeight());

		} else {
			currSiteStates = getNodeSequence(node);
//			// node is leaf or coalescent event
//			for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
//				int[] siteStates = ancestralTreeLikelihoodInput.get().get(i)
//						.getStatesForNode(ancestralTreeLikelihoodInput.get().get(0).treeInput.get(), node);
//
//				for (int j = 0; j < mapping.get(i).length; j++) {
//					currSiteStates[mapping.get(i)[j]] = siteStates[j];
//				}
//			}

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
			}else{
				if (startReadingFrame.get()!=null) {
					if (!getAminoAcid(startReadingFrame.get() - 1, currSiteStates).toString().equals("M"))
						System.err.println("AA sequence at root doesn't start with M + but with " + getAminoAcid(startReadingFrame.get() - 1, currSiteStates));

					if (endReadingFrame.get()==Integer.MAX_VALUE){
						// look for the first stop codon
						int i = startReadingFrame.get()-1;
						while (i<currSiteStates.length-2 && !isStopCodon(i, currSiteStates)){
							i=i+3;
						}

						// reading frame is inferred from the first stop codon
						endReadingFrame.set(i+1);
						System.err.println("the reading from is inferred to be from "+startReadingFrame.get()+" to "+ endReadingFrame.get() +
								", with " + (i+1) + " being the index of the stop codon");

					}
				}

			}
		}

		if (!node.isLeaf()) {
			if (node.getLeft()==null || node.getRight()==null)
				System.err.println("Node " + node.getHeight() + " has one child that is null");
			mapMuts(node.getLeft(), currSiteStates);
			mapMuts(node.getRight(), currSiteStates);
		}
	}

	protected int[] getNodeSequence(Node node){
		int[] currSiteStates = new int[totalLength];
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			int[] siteStates = ancestralTreeLikelihoodInput.get().get(i)
					.getStatesForNode(ancestralTreeLikelihoodInput.get().get(0).treeInput.get(), node);

			for (int j = 0; j < mapping.get(i).length; j++) {
				currSiteStates[mapping.get(i)[j]] = siteStates[j];
			}
		}
		return currSiteStates;
	}

	protected int[] getTranslatedNodeSequence(Node node){
		int[] currSiteStates = new int[totalLength];
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			int[] siteStates = ancestralTreeLikelihoodInput.get().get(i)
					.getStatesForNode(ancestralTreeLikelihoodInput.get().get(0).treeInput.get(), node);

			for (int j = 0; j < mapping.get(i).length; j++) {
				currSiteStates[mapping.get(i)[j]] = siteStates[j];
			}
		}

		int[] currTranslatedSites = new int[(endReadingFrame.get()-startReadingFrame.get())/3];

		for (int i = 0; i < currTranslatedSites.length; i++) {
			// get the index of the amino acid in the array aaaray
			currTranslatedSites[i] = aaList.indexOf(getAminoAcid(startReadingFrame.get()+i*3 - 1, currSiteStates));
		}

//		for (int j = 0; j < currTranslatedSites.length; j++){
//			System.out.print(aaList.get(currTranslatedSites[j]));
//		}
//		System.out.println();


		return currTranslatedSites;
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
	}

	private void initTable() {

		translationsTable = new HashMap<>();
		String[] dnaarray = dna.split(" ");
		String[] aaarray = aa.split(" ");
		for (int i = 0; i < dnaarray.length; i++)
			translationsTable.put(dnaarray[i], aaarray[i]);

		// make a list with the keys of translationTable
		aaList = Arrays.asList(aaarray).stream().distinct().collect(Collectors.toList());
		stopCodons = "TAG TAA TGA";

	}

	public String getMutation(int i, int p, int c) {
		return dataType.getCharacter(p) + "" + (i + 1) + ""	+ dataType.getCharacter(c);
	}

	public String getAAchange(int i, int[] parent, int[] child){
		if (i<startReadingFrame.get()-1 || i>endReadingFrame.get()-2){
			return "nan";
		}

		int codonPosition = (i-startReadingFrame.get()+1) % 3;

		int aaPosition = (i-startReadingFrame.get()+1) / 3;


		String codon_p = "";
		String codon_c = "";
		for (int j = 0; j < 3; j++) {
			codon_p = codon_p + dataType.getCharacter(parent[i+j-codonPosition]);
			codon_c = codon_c + dataType.getCharacter(child[i+j-codonPosition]);
		}
		return translationsTable.get(codon_p) + (aaPosition+1) + translationsTable.get(codon_c);
	}

	public String getAminoAcid(int i, int[] seq){
		if (i<startReadingFrame.get()-1 || i>endReadingFrame.get()-2){
			return "nan";
		}

		int codonPosition = (i-startReadingFrame.get()+1) % 3;
		String codon = "";
		for (int j = 0; j < 3; j++) {
			codon = codon + dataType.getCharacter(seq[i+j-codonPosition]);
		}
		return translationsTable.get(codon);
	}

	public boolean isStopCodon(int i, int[] seq){

		String codon = "";
		for (int j = 0; j < 3; j++) {
			codon = codon + dataType.getCharacter(seq[i+j]);
		}
		return stopCodons.contains(codon);
	}

	public String getBranchType(Node n){
		return mappedMascotInput.get().dynamics.getStringStateValue((int) n.getParent().getMetaData("location"))
				+ "->" + mappedMascotInput.get().dynamics.getStringStateValue((int) n.getMetaData("location"));
	}

}