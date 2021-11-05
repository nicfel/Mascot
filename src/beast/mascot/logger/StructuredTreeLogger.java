package beast.mascot.logger;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
//import org.jblas.DoubleMatrix;

import beast.core.Citation;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.StateNode;
import beast.core.parameter.BooleanParameter;

import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.coalescent.IntervalType;
import beast.mascot.distribution.Mascot;
import beast.mascot.ode.Euler2ndOrderTransitions;
import beast.mascot.ode.MascotODEUpDown;
import cern.colt.Arrays;


/**
 * @author Nicola Felix Mueller (nicola.felix.mueller@gmail.com)
 */
@Citation(	"Nicola F. Müller, David A. Rasmussen, Tanja Stadler (2018)\n"+
			"  MASCOT: parameter and state inference under the marginal\n"+
			"  structured coalescent approximation\n"+
			"  Bioinformatics, , bty406, https://doi.org/10.1093/bioinformatics/bty406")
public class StructuredTreeLogger extends Tree implements Loggable {
	
	
	public Input<Mascot> mascotInput = new Input<>("mascot", "Input of rates", Input.Validate.REQUIRED);
    
	public Input<Double> epsilonInput = new Input<>("epsilon", "step size for the RK4 integration",0.00001);
	public Input<Double> maxStepInput = new Input<>("maxStep", "step size for the RK4 integration", Double.POSITIVE_INFINITY);
	public Input<Double> stepSizeInput = new Input<>("stepSize", "step size for the RK4 integration");
	
	public Input<Boolean> useUpDown = new Input<>("upDown", "if up down algorithm is to use for the node state calculation", true);
    
    
    public Input<BranchRateModel.Base> clockModelInput = new Input<BranchRateModel.Base>("branchratemodel", "rate to be logged with branches of the tree");
    public Input<List<Function>> parameterInput = new Input<List<Function>>("metadata", "meta data to be logged with the tree nodes", new ArrayList<>());
    public Input<Boolean> maxStateInput = new Input<Boolean>("maxState", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<BooleanParameter> conditionalStateProbsInput = new Input<BooleanParameter>("conditionalStateProbs", "report branch lengths as substitutions (branch length times clock rate for the branch)");
    public Input<Boolean> substitutionsInput = new Input<Boolean>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<Integer> decimalPlacesInput = new Input<Integer>("dp", "the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)", -1);

    
    protected boolean someMetaDataNeedsLogging;
    protected boolean substitutions = false;
    protected boolean takeMax = true;
    protected boolean conditionals = true;
    
    protected boolean updown = true;

    protected DecimalFormat df;
    protected String type;
    
    protected int states;
    protected boolean[] used;
    protected boolean report;
    protected Mascot mascot;
    
	TreeInterface tree;

	
    @Override
    public void initAndValidate() {
    	mascot = mascotInput.get();
    	// RRB: correct?
    	tree = mascot.structuredTreeIntervalsInput.get().treeInput.get();
    	
        if (parameterInput.get().size() == 0 && clockModelInput.get() == null) {
        	someMetaDataNeedsLogging = false;
        	return;
            //throw new Exception("At least one of the metadata and branchratemodel inputs must be defined");
        }
    	someMetaDataNeedsLogging = true;
    	// without substitution model, reporting substitutions == reporting branch lengths 
        if (clockModelInput.get() != null) {
        	substitutions = substitutionsInput.get();
        }
       
        if (maxStateInput.get() != null){
        	takeMax = maxStateInput.get();
        	
        }

        int dp = decimalPlacesInput.get();

        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }             
        
        states = 0;    
    }

    @Override
    public void init(PrintStream out) {
    	mascot.structuredTreeIntervalsInput.get().treeInput.get().init(out);
    	states = mascot.dynamicsInput.get().getDimension();
   }

    public void log(int nSample, PrintStream out) {
    	log((long) nSample, out);
    }
    
    
    @Override
    public void log(long nSample, PrintStream out) {
    	states = mascotInput.get().dynamicsInput.get().getDimension();
    	
        // make sure we get the current version of the inputs
//        Tree tree = (Tree) mascot.treeIntervalsInput.get().treeInput.get().getCurrent();
        //calculate the state of each node
    	try {
			CalculateNodeStates();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    	used = new boolean[stateProbabilities.length];
    	report = false;
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode) {
        		metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        	}
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
        mascot.structuredTreeIntervalsInput.get().treeInput.get().getRoot().sort();
        root = mascot.structuredTreeIntervalsInput.get().treeInput.get().getRoot();
        out.print(toNewick(root, metadata, branchRateModel));
        out.print(";");
        
        for (int i = 0; i < used.length; i++)
        	if(!used[i])
        		System.err.println("not all nodes used");
        if (report)
        	System.err.println("error in node numbers");
    }
    Node root;

	/**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }

    String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
        if (maxStateInput.get() != null){
        	takeMax = maxStateInput.get();
        	
        }
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
        	if (leftID[node.getNr()-nrSamples] != node.getRight().getNr() && leftID[node.getNr()-nrSamples] != node.getLeft().getNr()){
        		report = true;
        		System.out.println("wrong nr of internal node: " 
        				+ leftID[node.getNr()-nrSamples] + " " + node.getLeft().getNr() + " " 
        				+ rightID[node.getNr()-nrSamples] + " " + node.getRight().getNr());
        		System.out.println(node.isRoot() + " " + node.getParent().isRoot() + " " + node.getLeft().getID() + " " + node.getRight().getID());
        		System.out.println(mascot.treeIntervalsInput.get().treeInput.get());
        		System.out.println(node.getTree());
        		System.exit(0);
        	}
        	if (!takeMax){	        
		        buf.append("[&");
		        
		        double[] stateProbs = getStateProb(node.getNr());		        
		        
		        for (int i = 0 ; i < states; i++)
		        	buf.append(String.format(Locale.US, "%s=%.3f,", mascot.dynamicsInput.get().getStringStateValue(i), stateProbs[i]));
		        
//		        buf.append(String.format("%.3f", stateProbs.get(states-1)));
//		        buf.append("}");
	        
		        buf.append("max=");
		        buf.append(String.format("%s", 
		        		mascot.dynamicsInput.get().getStringStateValue(whichMax(stateProbs))));
		        
//		        if (node.getLeft().isLeaf()){
//		        	buf.append(String.format("left=%d,", (int) mascot.dynamicsInput.get().typeTraitInput.get().getValue(node.getLeft().getID() )) );
//		        }else{
//		        	buf.append("left=Nan,");
//		        }
//		        if (node.getRight().isLeaf()){
//		        	buf.append(String.format("right=%d", (int) mascot.dynamicsInput.get().typeTraitInput.get().getValue(node.getRight().getID() )) );
//		        }else{
//		        	buf.append("right=Nan");
//		        }
		        	
		        
		        if (branchRateModel != null) {
		            buf.append(",rate=");
	                appendDouble(buf, branchRateModel.getRateForBranch(node));
		        }
		        buf.append(']');
        	}else{
		        buf.append("[&max" + type + "=");
		        double[] stateProbs = getStateProb(node.getNr());		        

		        buf.append(String.format("%d", whichMax(stateProbs) ));
		        
		        if (branchRateModel != null) {
		            buf.append(",rate=");
	                appendDouble(buf, branchRateModel.getRateForBranch(node));
		        }

		        buf.append(']');        		
        	}
        }else{
			String sampleID = node.getID();
			String[] splits = sampleID.split("_");
			int sampleState;
			
			if(mascot.dynamicsInput.get().typeTraitInput.get()!=null){				
				sampleState = mascot.dynamicsInput.get().getValue(node.getID());
			}			
			
			else{
				sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
			}
			if (!takeMax){
    	        
		        buf.append("[&");
	
		        for (int i = 0 ; i < states; i++){
		        	if (sampleState != i) buf.append(String.format("%s=0,", mascot.dynamicsInput.get().getStringStateValue(i)));
		        	if (sampleState == i) buf.append(String.format("%s=1,", mascot.dynamicsInput.get().getStringStateValue(i)));
		        }		        
		        buf.append("max=");

		        buf.append(String.format("%s", 
		        		mascot.dynamicsInput.get().getStringStateValue(sampleState)) );
		        buf.append(']');
        	}else{
		        buf.append("[&max" + type + "=");

		        buf.append(String.format("%d", sampleState ));
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
		mascot.structuredTreeIntervalsInput.get().treeInput.get().close(out);
    }
	
	//===================================================
	//===================================================
	// Calculate the state of all nodes using the up-down
	// algorithm
	//===================================================
	//===================================================
	public int samples;
	public int nrSamples;
	public double[][] stateProbabilities;
	public double[][] stateProbabilitiesDown;
	public double[][][] TransitionProbabilities;	  
	public int[] leftID;
	public int[] rightID;
    
    public int nrLineages;   

    // current rates         
    private double[] migrationRates;
    private int[] indicators;
    protected double[] coalescentRates; 	

    
    // Set up for lineage state probabilities
    protected ArrayList<Integer> activeLineages;
	private double[] linProbs;
	private double[] transitionProbs;
	
	// maximum integration error tolerance
    private double maxTolerance = 1e-5;            
    private boolean recalculateLogP;   

	
    public void CalculateNodeStates() throws Exception{  
    	// newly calculate tree intervals
    	mascot.structuredTreeIntervalsInput.get().calculateIntervals();
    	// correctly calculate the daughter nodes at coalescent intervals in the case of
    	// bifurcation or in case two nodes are at the same height
    	mascot.structuredTreeIntervalsInput.get().swap();    	
 
    	leftID = new int[mascot.structuredTreeIntervalsInput.get().getSampleCount()];
    	rightID = new int[mascot.structuredTreeIntervalsInput.get().getSampleCount()];

    	
    	stateProbabilities = new double[mascot.structuredTreeIntervalsInput.get().getSampleCount()][];
    	stateProbabilitiesDown = new double[mascot.structuredTreeIntervalsInput.get().getSampleCount()][];
    	TransitionProbabilities = new double[mascot.structuredTreeIntervalsInput.get().getSampleCount()*2][][];        
        nrSamples = mascot.structuredTreeIntervalsInput.get().getSampleCount() + 1;        

    	
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages = new ArrayList<Integer>(); 

        nrLineages = 0;
        linProbs = new double[0];// initialize the tree and rates interval counter
        transitionProbs = new double[0];// initialize the tree and rates interval counter
        
        int treeInterval = 0, ratesInterval = 0;        
        double nextEventTime = 0.0;
		coalescentRates = mascot.dynamicsInput.get().getCoalescentRate(ratesInterval);  
        migrationRates = mascot.dynamicsInput.get().getBackwardsMigration(ratesInterval);
		indicators = mascot.dynamicsInput.get().getIndicators(ratesInterval);  
        // Time to the next rate shift or event on the tree
        double nextTreeEvent = mascot.structuredTreeIntervalsInput.get().getInterval(treeInterval);
        double nextRateShift = mascot.dynamicsInput.get().getInterval(ratesInterval);

             

        int currTreeInterval = 0; 													// what tree interval are we in?
        // Calculate the likelihood
        do {       
        	nextEventTime = Math.min(nextTreeEvent, nextRateShift); 	 
        	
        	if (nextEventTime > 0) {													// if true, calculate the interval contribution        		
                if(recalculateLogP){
    				System.err.println("ode calculation stuck, reducing tolerance, new tolerance= " + maxTolerance);
    				maxTolerance *=0.1;
                	CalculateNodeStates();
                	return;
                }
                if(stepSizeInput.get()!=null){
    	        	double[] probs_for_ode = new double[linProbs.length + transitionProbs.length];   
    	        	double[] oldLinProbs = new double[linProbs.length + transitionProbs.length]; 
    	        	
                    for (int i = 0; i<linProbs.length; i++)
                    	oldLinProbs[i] = linProbs[i];  
                    for (int i = linProbs.length; i < transitionProbs.length; i++)
                    	oldLinProbs[i] = transitionProbs[i-linProbs.length];  

                	
	                FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(stepSizeInput.get());	                
//	                integrator.setMaxEvaluations((int) 1e5);  // set the maximal number of evaluations              
	                FirstOrderDifferentialEquations ode = new MascotODEUpDown(migrationRates, coalescentRates, nrLineages , coalescentRates.length);

	                // integrate	        
	                try {
	                	integrator.integrate(ode, 0, oldLinProbs, nextEventTime, probs_for_ode);
	                }catch(Exception e){
	                	System.out.println(e);
	                	System.out.println("expection");
	                	System.exit(0);
	                	recalculateLogP = true;    				
	                }        		       	
		           
	                for (int i = 0; i<linProbs.length; i++)
	            		linProbs[i] = probs_for_ode[i];  
	                for (int i = linProbs.length; i < transitionProbs.length; i++)
	            		transitionProbs[i-linProbs.length] = probs_for_ode[i];  
	        	}else {
		        	double[] linProbs_tmp = new double[linProbs.length + transitionProbs.length]; 
		        	double[] linProbs_tmpdt = new double[linProbs.length + transitionProbs.length]; 
		        	double[] linProbs_tmpddt = new double[linProbs.length + transitionProbs.length]; 
		        	double[] linProbs_tmpdddt = new double[linProbs.length + transitionProbs.length]; 

                    for (int i = 0; i<linProbs.length; i++)
                    	linProbs_tmp[i] = linProbs[i];  
                    
                    for (int i = linProbs.length; i < (transitionProbs.length+linProbs.length); i++)
                    	linProbs_tmp[i] = transitionProbs[i-linProbs.length];   		
    	        	
                    Euler2ndOrderTransitions euler;
	        		if (mascot.dynamicsInput.get().hasIndicators)
	        			euler = new Euler2ndOrderTransitions(migrationRates, indicators, coalescentRates, nrLineages , coalescentRates.length, epsilonInput.get(), maxStepInput.get());
	        		else
	        			euler = new Euler2ndOrderTransitions(migrationRates, coalescentRates, nrLineages , coalescentRates.length, epsilonInput.get(), maxStepInput.get());
		        	
		        	
		        	linProbs[linProbs.length-1] = 0;
		        	euler.calculateValues(nextEventTime, linProbs_tmp, linProbs_tmpdt, linProbs_tmpddt, linProbs_tmpdddt);		        	
	        		
	                for (int i = 0; i<linProbs.length; i++)
	            		linProbs[i] = linProbs_tmp[i];  
	                for (int i = linProbs.length; i < linProbs_tmp.length; i++)
	            		transitionProbs[i-linProbs.length] = linProbs_tmp[i];  

	        	}
			}
        	
        	if (nextTreeEvent <= nextRateShift){
 	        	if (mascot.structuredTreeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.COALESCENT) {
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		normalizeLineages();									// normalize all lineages before event		
        			coalesce(treeInterval);	  				// calculate the likelihood of the coalescent event
	        	}
 	       		
 	       		if (mascot.structuredTreeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	
 	       			nrLineages++;													// sampling event increases the number of lineages by one
 	       			if (linProbs.length > 0)
 	       				normalizeLineages();								// normalize all lineages before event
 	       			sample(treeInterval);							// calculate the likelihood of the sampling event if sampling rate is given
	       		}	
 	       		
 	       		treeInterval++;
        		nextRateShift -= nextTreeEvent;   
        		try{
        			nextTreeEvent = mascot.structuredTreeIntervalsInput.get().getInterval(treeInterval);
        		}catch(Exception e){
        			break;
        		}
        	}else{
        		ratesInterval++;
        		coalescentRates = mascot.dynamicsInput.get().getCoalescentRate(ratesInterval);  
                migrationRates = mascot.dynamicsInput.get().getBackwardsMigration(ratesInterval);
        		indicators = mascot.dynamicsInput.get().getIndicators(ratesInterval);  
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = mascot.dynamicsInput.get().getInterval(ratesInterval);
        	}
        	
        }while(nextTreeEvent <= Double.POSITIVE_INFINITY);
        currTreeInterval = mascot.structuredTreeIntervalsInput.get().getIntervalCount()-1;
        
        do{
		  	if (mascot.structuredTreeIntervalsInput.get().getIntervalType(currTreeInterval) == IntervalType.COALESCENT) {
		  		coalesceDown(currTreeInterval);									// Set parent lineage state probs and remove children
		   	}       	
		  	currTreeInterval--;
        }while(currTreeInterval>=0);
//        System.exit(0);
    }
	

    private double normalizeLineages(){
    	if (linProbs==null)
    		return 0.0;
    	
    	
    	double interval = 0.0;
    	for (int i = 0; i < linProbs.length/states; i++){
    		double lineProbs = 0.0;
    		for (int j = 0; j < states; j++)
    			if (linProbs[i*states+j]>=0.0){
    				lineProbs += linProbs[i*states+j];
    			}else{
    				// try recalculation after lowering the tolerance
    				recalculateLogP = true;
    				return Math.log(1.0);
    			}
    		for (int j = 0; j < states; j++){
    			linProbs[i*states+j] = linProbs[i*states+j]/lineProbs;
    		}    		
    		interval +=lineProbs;
    	}	
    	
		// return mean P_t(T)
		return Math.log(interval/(linProbs.length/states));

    }
    
    private void sample(int currTreeInterval) {
		int incomingLines = mascot.structuredTreeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		// calculate the new length of the arrays for the transition and lineage states
		int newLengthLineages = linProbs.length + 1*states;
		int newLengthTransitions = transitionProbs.length + 1*states*states;
		
		double[] linProbsNew = new double[newLengthLineages];		
		double[] transitionProbsNew = new double[newLengthTransitions];		
		
		for (int i = 0; i < linProbs.length; i++)
			linProbsNew[i] = linProbs[i];	
		
		for (int i = 0; i < transitionProbs.length; i++)
			transitionProbsNew[i] = transitionProbs[i];		
		
		int currPositionLineages = linProbs.length;
		int currPositionTransitions = transitionProbs.length;
		/*
		 * If there is no trait given as Input, the model will simply assume that
		 * the last value of the taxon name, the last value after a _, is an integer
		 * that gives the type of that taxon
		 */
		if (mascot.dynamicsInput.get().typeTraitInput.get()!=null){
			int  l = incomingLines; {
				activeLineages.add(l);//.getNr());
				int sampleState = (int) mascot.dynamicsInput.get().getValue(tree.getNode(l).getID());
				for (int i = 0; i< states; i++){
					if (i == sampleState){
						linProbsNew[currPositionLineages] = 1.0;currPositionLineages++;
					}
					else{
						linProbsNew[currPositionLineages] = 0.0;currPositionLineages++;
					}
				}
				// add the initial transition probabilities (diagonal matrix)
				for (int s = 0; s < states; s++){
					for (int i = 0; i < states; i++){
						if (i == s){
							transitionProbsNew[currPositionTransitions] = 1.0;
							currPositionTransitions++;
						}else{
							transitionProbsNew[currPositionTransitions] = 0.0;
							currPositionTransitions++;
						}
					}
				}
			}				
		}else{
			int l = incomingLines; {
				activeLineages.add(l);//.getNr());
				String sampleID = tree.getNode(l).getID();
				int sampleState = 0;
				if (states > 1){				
					String[] splits = sampleID.split("_");
					sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
				}
				for (int i = 0; i< states; i++){
					if (i == sampleState){
						linProbsNew[currPositionLineages] = 1.0;currPositionLineages++;
					}
					else{
						linProbsNew[currPositionLineages] = 0.0;currPositionLineages++;
					}
				}
				// add the initial transition probabilities (diagonal matrix)
				for (int s = 0; s < states; s++){
					for (int i = 0; i < states; i++){
						if (i == s){
							transitionProbsNew[currPositionTransitions] = 1.0;
							currPositionTransitions++;
						}else{
							transitionProbsNew[currPositionTransitions] = 0.0;
							currPositionTransitions++;
						}
					}
				}
			}	
		}
		linProbs = linProbsNew;
		transitionProbs = transitionProbsNew;
    }
    
    private void coalesce(int currTreeInterval) {
    	// normalize the transition probabilities
    	for (int i = 0; i < nrLineages*states; i++){
    		double lineProbs = 0.0;
    		for (int j = 0; j < states; j++){
    			if (transitionProbs.length>=0.0){
    				lineProbs += transitionProbs[i*states+j];
    			}else{
    				System.err.println("transition probability smaller than 0 (or NaN before normalizing");	    				
    				System.exit(0);
    			}
    		}
    		for (int j = 0; j < states; j++)
    			transitionProbs[i*states+j] = transitionProbs[states*i+j]/lineProbs; 	    		
    	}

    	

//		List<Integer> coalLines = mascot.treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
//    	if (coalLines.size() > 2) {
//			System.err.println("Unsupported coalescent at non-binary node");
//			System.exit(0);
//		}
//    	if (coalLines.size() < 2) {
//    		System.out.println();
//    		System.out.println("WARNING: Less than two lineages found at coalescent event!");
//    		System.exit(0);
//		}
		
    	int [] coalLines = new int[] {
    			mascot.structuredTreeIntervalsInput.get().getLineagesRemoved(currTreeInterval,0),
    			mascot.structuredTreeIntervalsInput.get().getLineagesRemoved(currTreeInterval,1)
    	};
    	final int daughterIndex1 = activeLineages.indexOf(coalLines[0]);//.getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines[1]);//.getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println("daughter lineages at coalescent event not found");
    		System.exit(0);
		}
		double[] lambda = new double[states];
		double lambdaSum = 0;

		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */
        for (int k = 0; k < states; k++) { 
        	Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];			
			if (!Double.isNaN(pairCoalRate)){
				lambda[k] =  pairCoalRate;
				lambdaSum += pairCoalRate;
			}else{
//	    		System.exit(0);
			}
        }
        
        activeLineages.add(tree.getNode(coalLines[0]).getParent().getNr());
        
        
        // get the node state probabilities
		double[] pVec = new double[states];
		for (int i = 0; i < pVec.length; i++)
			pVec[i] = lambda[i]/lambdaSum;
		
		// save the node states conditioned on the subtree
		stateProbabilities[tree.getNode(coalLines[0]).getParent().getNr() - nrSamples] = pVec;

		
		// get the transition probabilities of daughter lineage 1
		double[][] tP1 = new double[states][states];
		
		for (int i = 0; i< states; i++){
			for (int j = 0; j< states; j++){
				tP1[i][j] = transitionProbs[daughterIndex1*states*states+i*states+j];
//				tP1.put(i, j, transitionProbs[daughterIndex1*states*states+i*states+j]);
			}
		}
//		tP1.print();
//		if (!coalLines.get(0).isLeaf()){
//			DoubleMatrix start = stateProbabilities[coalLines.get(0).getNr() - nrSamples];
////			start.print();
//			start.transpose().mmul(tP1).div(start.transpose().mmul(tP1).sum()).print();
//			System.out.print("[");
//			for (int i = 0; i < (states-1); i++)
//				System.out.print(String.format("%.6f", linProbs[daughterIndex1*states+i]) + ", ");
//			System.out.print(String.format("%.6f", linProbs[daughterIndex1*states+states-1]) + "]\n");
//			System.out.println();
////			System.exit(0);
//		}
		
		// get the transition probabilities of daughter lineage 2
		double[][] tP2 = new double[states][states];
		for (int i = 0; i< states; i++){
			for (int j = 0; j< states; j++){
				tP2[i][j] = transitionProbs[daughterIndex2*states*states+i*states+j];
//				tP2.put(i, j, transitionProbs[daughterIndex2*states*states+i*states+j]);
			}
		}	
		
		double[] linProbsNew  = new double[linProbs.length - states];

		int linCount = 0;		
		// add all lineages execpt the daughter lineage to the new p array
		for (int i = 0; i <= nrLineages; i++){
			if (i != daughterIndex1 && i != daughterIndex2){
				for (int j = 0; j < states; j++){
					linProbsNew[linCount*states + j] = linProbs[i*states + j];
				}
				linCount++;
			}
		}
		// add the parent lineage
		for (int j = 0; j < states; j++){
			linProbsNew[linCount*states + j] = pVec[j];
		}
		// set p to pnew
		linProbs = linProbsNew;	

		
		double[] transitionProbsNew  = new double[transitionProbs.length - states*states];
		
		// add initial transition probabilities for the parent lineage
		linCount = 0;
		for (int i = 0; i <= nrLineages; i++){
			if (i != daughterIndex1 && i != daughterIndex2){
				for (int j = 0; j < states; j++)
					for (int k = 0; k < states; k++)
						transitionProbsNew[linCount*states*states+j*states+k] 
								= transitionProbs[i*states*states+j*states+k];
				linCount++;
			}
		}

		for (int j = 0; j < states; j++)
			for (int k = 0; k < states; k++)
				if (j==k)
					transitionProbsNew[linCount*states*states+j*states+k] = 1.0;
				else
					transitionProbsNew[linCount*states*states+j*states+k] = 0.0;

		// set the transition probs
		transitionProbs = transitionProbsNew;
		
		// save the transition probabilities of each of the two daughter lineages
		TransitionProbabilities[coalLines[0]/*.getNr()*/] = tP1;
		TransitionProbabilities[coalLines[1]/*.getNr()*/] = tP2;
		
		//Remove daughter lineages from the line state probs and the transition probs
		if (daughterIndex1>daughterIndex2){
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex1);
			activeLineages.remove(daughterIndex2);			
		}else{
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex2);
			activeLineages.remove(daughterIndex1);
		}
		
		if(tree.getNode(coalLines[0]).getParent().getNr() != mascot.structuredTreeIntervalsInput.get().treeInput.get().getNode(coalLines[1]/*.getNr()*/).getParent().getNr())
			System.err.println("wrong daughter parent");
		if(tree.getNode(coalLines[1]).getParent().getNr() != mascot.structuredTreeIntervalsInput.get().treeInput.get().getNode(coalLines[0]/*.getNr()*/).getParent().getNr())
			System.err.println("wrong daughter parent");
		if(tree.getNode(coalLines[1]).getParent().getNr() != tree.getNode(coalLines[0]).getParent().getNr())
			System.err.println("coalescent nodes don't have the same parent");
		
		leftID[tree.getNode(coalLines[0]).getParent().getNr() - nrSamples] = tree.getNode(coalLines[0]).getParent().getLeft().getNr();
		rightID[tree.getNode(coalLines[0]).getParent().getNr() - nrSamples] = tree.getNode(coalLines[0]).getParent().getRight().getNr();
    }

    private void coalesceDown(int currTreeInterval) {
		int parentLines = mascot.structuredTreeIntervalsInput.get().getLineagesAdded(currTreeInterval);
//		if (parentLines.size()!=1){
//			System.err.println("to many lineages, while coalescening down");
//			System.exit(0);
//		}
		Node parentNode = tree.getNode(parentLines);
		
		if (!parentNode.isRoot()){
			double[] start = stateProbabilities[parentNode.getNr() - nrSamples];
			double[] end = stateProbabilitiesDown[parentNode.getParent().getNr() - nrSamples];
			double[][] flow = TransitionProbabilities[parentNode.getNr()];
			double[] otherSideInfo = new double[states];
			for (int a = 0; a < states; a++) {
				double sum = 0;
				for (int b = 0; b < states; b++) { 
					sum += start[b] * flow[b][a];
				}
				otherSideInfo[a] = end[a]/sum;
				if (Double.isNaN(otherSideInfo[a]))
					otherSideInfo[a] = 0;
//				DoubleMatrix otherSideInfo = end.div(start.transpose().mmul(flow));
			}
//			// get rid of NaN from division by 0
//			for (int i = 0; i < otherSideInfo.length; i++)
//				if (Double.isNaN(otherSideInfo.get(i)))
//					otherSideInfo.put(i, 0.0);

			double[] conditional = new double[states];
			double condsum = 0;
			for (int a = 0; a < states; a++) {
				double sum = 0;
				for (int b = 0; b < states; b++) { 
					sum += flow[a][b] * otherSideInfo[b];
				}
				conditional[a] = sum *start[a];
				condsum += conditional[a];
			}
			for (int a = 0; a < states; a++)
				conditional[a] /= condsum;


//			DoubleMatrix conditional = flow.mmul(otherSideInfo);
			
//			conditional = conditional.mul(start);
			stateProbabilitiesDown[parentNode.getNr() - nrSamples] = conditional;
//			if (!(conditional.get(0) >= 0.0 && conditional.get(0)<=1.0))
//				conditional.print();
		}else{
			stateProbabilitiesDown[parentNode.getNr() - nrSamples] = stateProbabilities[parentNode.getNr() - nrSamples];
    	}
	}

    
	public double[] getStateProb(int nr) {
		if(useUpDown.get()){			
			used[nr - nrSamples] = true;
			return stateProbabilitiesDown[nr - nrSamples] ;
		}else{
			used[nr - nrSamples] = true;
			return stateProbabilities[nr - nrSamples] ;
		}
	}
	
	public double[] getStateProbOnly(int nr) {
		if(useUpDown.get()){			
			return stateProbabilitiesDown[nr - nrSamples] ;
		}else{
			return stateProbabilities[nr - nrSamples] ;
		}
	}
	
	public int whichMax(double[] stateProbs) {
		double max_val = -1;
		int max_ind = 1;
		for (int i = 0; i < stateProbs.length;i++) {
			if (stateProbs[i]>max_val) {
				max_val = stateProbs[i];
				max_ind = i;
			}
		}
		return max_ind;
	}

	
    public void calcForTest() {
    	states = mascotInput.get().dynamicsInput.get().getDimension();
    	
    	
    	try {
			CalculateNodeStates();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

	
}
