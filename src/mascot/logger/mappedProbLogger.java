package mascot.logger;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.BooleanParameter;
import mascot.distribution.MappedMascot;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Nicola Felix Mueller
 */
@Description("logs the state of the root, i.e. the probability the root being in "+ 
			"any of m states based on the MASCOT density class")
public class mappedProbLogger extends CalculationNode implements Loggable {
	public Input<StructuredTreeLogger> stlInput = new Input<StructuredTreeLogger>(
			"structuredTreeLogger", "stl",Input.Validate.REQUIRED);
	public Input<MappedMascot> mappedInput = new Input<MappedMascot>(
			"mappedTreeLogger", "stl",Input.Validate.REQUIRED);
	
	
	public Input<BranchRateModel.Base> clockModelInput = new Input<BranchRateModel.Base>("branchratemodel",
			"rate to be logged with branches of the tree");
	public Input<List<Function>> parameterInput = new Input<List<Function>>("metadata",
			"meta data to be logged with the tree nodes", new ArrayList<>());
	public Input<Boolean> maxStateInput = new Input<Boolean>("maxState",
			"report branch lengths as substitutions (branch length times clock rate for the branch)", false);
	public Input<BooleanParameter> conditionalStateProbsInput = new Input<BooleanParameter>("conditionalStateProbs",
			"report branch lengths as substitutions (branch length times clock rate for the branch)");
	public Input<Boolean> substitutionsInput = new Input<Boolean>("substitutions",
			"report branch lengths as substitutions (branch length times clock rate for the branch)", false);
	public Input<Integer> decimalPlacesInput = new Input<Integer>("dp",
			"the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)",
			-1);


	
	private int states;
	private Node root;
	protected DecimalFormat df;
	protected boolean someMetaDataNeedsLogging;
	protected boolean substitutions = false;


	@Override
	public void initAndValidate() {
		if (parameterInput.get().size() == 0 && clockModelInput.get() == null) {
			someMetaDataNeedsLogging = false;
			return;
			// throw new Exception("At least one of the metadata and branchratemodel inputs
			// must be defined");
		}
		someMetaDataNeedsLogging = true;
		// without substitution model, reporting substitutions == reporting branch
		// lengths
		if (clockModelInput.get() != null) {
			substitutions = substitutionsInput.get();
		}

		int dp = decimalPlacesInput.get();

		if (dp < 0) {
			df = null;
		} else {
			// just new DecimalFormat("#.######") (with dp time '#' after the decimal)
			df = new DecimalFormat("#." + new String(new char[dp]).replace('\0', '#'));
			df.setRoundingMode(RoundingMode.HALF_UP);
		}
	}


	
	@Override
	public void init(PrintStream out) {
		stlInput.get().init(out);
	}

	@Override
	public void log(long nSample, PrintStream out) {
		out.print("tree STATE_" + nSample + " = ");
		
		mappedInput.get().calculateLogP();
//		mappedInput.get().log(0, System.out);

		stlInput.get().calcForTest();
		
		root = mappedInput.get().getRoot();
		out.print(toNewick(root, null, null, stlInput.get()));
		out.print(";");		
	}
	
	private void appendDouble(StringBuffer buf, double d) {
		if (df == null) {
			buf.append(d);
		} else {
			buf.append(df.format(d));
		}
	}

	
	String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel, StructuredTreeLogger stl) {
		StringBuffer buf = new StringBuffer();
		if (node.getLeft() != null) {
			buf.append("(");
			buf.append(toNewick(node.getLeft(), metadataList, branchRateModel, stlInput.get()));
			if (node.getRight() != null) {
				buf.append(',');
				buf.append(toNewick(node.getRight(), metadataList, branchRateModel, stlInput.get()));
			}
			buf.append(")");
		} else {
			buf.append(node.getNr() + 1);
		}
		if (!node.isLeaf()) {
//			System.out.println(node.getMetaData("location"));;
//			System.out.println(node.getHeight());
//			System.out.println(lalala);
			
			buf.append("[&");
			buf.append(mappedInput.get().dynamics.typeTraitInput.getName() + "=" + mappedInput.get().dynamics.getStringStateValue((int) node.getMetaData("location")));
			
		
			if (branchRateModel != null) {
				buf.append(",rate=");
				appendDouble(buf, branchRateModel.getRateForBranch(node));
			}
			
			if (node.getChildCount()==2) {
//				System.out.println(node.getNr());
//				System.out.println(node.getHeight());
//				stl.log(0, System.out);
				double[] probs = stl.getStateProbOnly(node.getNr());
				for (int j = 0; j < probs.length; j++) {
					buf.append("," + mappedInput.get().dynamics.getStringStateValue(j) + "=");
					buf.append(probs[j]);					
				}
			}
			
			buf.append(']');


		} else {
			String sampleID = node.getID();
			String[] splits = sampleID.split("_");
			int sampleState;

			if (mappedInput.get().dynamics.typeTraitInput.get() != null) {
				sampleState = mappedInput.get().dynamics.getValue(node.getID());
			}

			else {
				sampleState = Integer.parseInt(splits[splits.length - 1]); // samples states (or priors) should
																			// eventually be specified in the XML
			}
			
			if ( node.getMetaData("location")!=null) {
			buf.append("[&");
			buf.append(mappedInput.get().dynamics.typeTraitInput.getName() + "=" + mappedInput.get().dynamics.getStringStateValue((int) node.getMetaData("location")));
			buf.append(']');
			}

		}

		buf.append(":");
		if (substitutions) {
			appendDouble(buf, node.getLength() * branchRateModel.getRateForBranch(node));
		} else {
			appendDouble(buf, node.getLength());
		}
		return buf.toString();
	}



	@Override
	public void close(PrintStream out) {	
		stlInput.get().close(out);

	}

	

	
}
