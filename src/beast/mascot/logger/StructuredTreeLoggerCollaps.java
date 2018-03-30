package beast.mascot.logger;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.jblas.DoubleMatrix;

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
import beast.mascot.distribution.MascotCollaps;
import beast.mascot.ode.Euler2ndOrder;
import beast.mascot.ode.Euler2ndOrderTransitions;
import beast.mascot.ode.MascotODE;
import beast.mascot.ode.MascotODEUpDown;


/**
 * @author Nicola Felix Mueller (nicola.felix.mueller@gmail.com)
 */
public class StructuredTreeLoggerCollaps extends Tree implements Loggable {
	
	
	public Input<MascotCollaps> mascotInput = new Input<>("mascot", "Input of rates", Input.Validate.REQUIRED);
    
	public Input<Double> epsilonInput = new Input<>("epsilon", "step size for the RK4 integration",0.00001);
	public Input<Double> maxStepInput = new Input<>("maxStep", "step size for the RK4 integration", Double.POSITIVE_INFINITY);
	public Input<Double> stepSizeInput = new Input<>("stepSize", "step size for the RK4 integration");
	    
    
    public Input<BranchRateModel.Base> clockModelInput = new Input<BranchRateModel.Base>("branchratemodel", "rate to be logged with branches of the tree");
    public Input<List<Function>> parameterInput = new Input<List<Function>>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());
    public Input<Boolean> maxStateInput = new Input<Boolean>("maxState", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<BooleanParameter> conditionalStateProbsInput = new Input<BooleanParameter>("conditionalStateProbs", "report branch lengths as substitutions (branch length times clock rate for the branch)");
    public Input<Boolean> substitutionsInput = new Input<Boolean>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<Integer> decimalPlacesInput = new Input<Integer>("dp", "the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)", -1);

    
    boolean someMetaDataNeedsLogging;
    boolean substitutions = false;
    boolean takeMax = true;
    boolean conditionals = true;
    
    boolean updown = false;

    private DecimalFormat df;
    private String type;
    
    private int states;
    boolean[] used;
	boolean report;
    

	TreeInterface tree;
	
    @Override
    public void initAndValidate() { 
    	// RRB: correct?
    	tree = this;
    	
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
    	mascotInput.get().treeIntervalsInput.get().treeInput.get().init(out);
    	states = mascotInput.get().dynamicsInput.get().getDimension();
   }

    @Override
    public void log(int nSample, PrintStream out) {
    	states = mascotInput.get().dynamicsInput.get().getDimension();
    	
        // make sure we get the current version of the inputs
//        Tree tree = (Tree) mascotInput.get().treeIntervalsInput.get().treeInput.get().getCurrent();
        //calculate the state of each node
    	calculateNodeStates();
    	
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
        mascotInput.get().treeIntervalsInput.get().treeInput.get().getRoot().sort();
        root = mascotInput.get().treeIntervalsInput.get().treeInput.get().getRoot();
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
        		System.out.println(mascotInput.get().treeIntervalsInput.get().treeInput.get());
        		System.out.println(node.getTree());
        		System.exit(0);
        	}
        	if (!takeMax){	        
		        buf.append("[&");
		        
		        DoubleMatrix stateProbs = new DoubleMatrix();
		        
	        	stateProbs = getStateProb(node.getNr());		        
		        if (node.getHeight() < mascotInput.get().collapsTimeInput.get().getValue()){
			        for (int i = 0 ; i < states; i++)
			        	buf.append(String.format("%s=%.3f,", mascotInput.get().dynamicsInput.get().getStringStateValue(i), stateProbs.get(i)));
		        	buf.append(String.format("ancestral=0.0,"));
			        buf.append("max=");
			        buf.append(String.format("%s", 
			        		mascotInput.get().dynamicsInput.get().getStringStateValue(stateProbs.argmax())));
		        }else{
			        for (int i = 0 ; i < states; i++)
		        		buf.append(String.format("%s=0.0,", mascotInput.get().dynamicsInput.get().getStringStateValue(i)));		        	
		        	buf.append(String.format("ancestral=1.0,"));
			        buf.append("max=");
			        buf.append(String.format("ancestral"));
		        }
		        
		        
//		        buf.append(String.format("%.3f", stateProbs.get(states-1)));
//		        buf.append("}");
	        
//		        buf.append("max=");
//		        buf.append(String.format("%s", 
//		        		mascotInput.get().dynamicsInput.get().getStringStateValue(stateProbs.argmax())));
		        
//		        if (node.getLeft().isLeaf()){
//		        	buf.append(String.format("left=%d,", (int) mascotInput.get().dynamicsInput.get().typeTraitInput.get().getValue(node.getLeft().getID() )) );
//		        }else{
//		        	buf.append("left=Nan,");
//		        }
//		        if (node.getRight().isLeaf()){
//		        	buf.append(String.format("right=%d", (int) mascotInput.get().dynamicsInput.get().typeTraitInput.get().getValue(node.getRight().getID() )) );
//		        }else{
//		        	buf.append("right=Nan");
//		        }
		        	
		        
		        
		        buf.append(']');
        	}else{
		        buf.append("[&max" + type + "=");
		        DoubleMatrix stateProbs = new DoubleMatrix();
		        
	        	stateProbs = getStateProb(node.getNr());

		        buf.append(String.format("%d", stateProbs.argmax() ));
		        buf.append(']');        		
        	}
        }else{
			String sampleID = node.getID();
			String[] splits = sampleID.split("_");
			int sampleState;
			
			if(mascotInput.get().dynamicsInput.get().typeTraitInput.get()!=null){				
				sampleState = mascotInput.get().dynamicsInput.get().getValue(node.getID());
			}			
			
			else{
				sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
			}
			if (!takeMax){
    	        
		        buf.append("[&");
	
		        for (int i = 0 ; i < states; i++){
		        	if (sampleState != i) buf.append(String.format("%s=0,", mascotInput.get().dynamicsInput.get().getStringStateValue(i)));
		        	if (sampleState == i) buf.append(String.format("%s=1,", mascotInput.get().dynamicsInput.get().getStringStateValue(i)));
		        }		        
	        	buf.append(String.format("ancestral=0.0,"));
		        buf.append("max=");

		        buf.append(String.format("%s", 
		        		mascotInput.get().dynamicsInput.get().getStringStateValue(sampleState)) );
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
		mascotInput.get().treeIntervalsInput.get().treeInput.get().close(out);
    }
	
	//===================================================
	//===================================================
	// Calculate the state of all nodes using the up-down
	// algorithm
	//===================================================
	//===================================================
	public int samples;
	public int nrSamples;
	public DoubleMatrix[] stateProbabilities;
	public DoubleMatrix[] stateProbabilitiesDown;
	public DoubleMatrix[] TransitionProbabilities;	  
	public int[] leftID;
	public int[] rightID;
    
    private int nrLineages;   

    // current rates         
    private double[] migrationRates;
    private int[] indicators;
    private double[] coalescentRates; 	

    
    // Set up for lineage state probabilities
    private ArrayList<Integer> activeLineages;
	private double[] linProbs;
	private double[] transitionProbs;
	
	// maximum integration error tolerance
    private double maxTolerance = 1e-5;            
    private boolean recalculateLogP;   

	
    public double calculateNodeStates() {
    	stateProbabilities = new DoubleMatrix[mascotInput.get().treeIntervalsInput.get().getSampleCount()];
        nrSamples = mascotInput.get().treeIntervalsInput.get().getSampleCount() + 1;    
    	leftID = new int[mascotInput.get().treeIntervalsInput.get().getSampleCount()];
    	rightID = new int[mascotInput.get().treeIntervalsInput.get().getSampleCount()];
    	// newly calculate tree intervals
    	mascotInput.get().treeIntervalsInput.get().calculateIntervals();
    	// correctly calculate the daughter nodes at coalescent intervals in the case of
    	// bifurcation or in case two nodes are at the same height
    	mascotInput.get().treeIntervalsInput.get().swap();    	
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages = new ArrayList<Integer>(); 
        double logP = 0;
        nrLineages = 0;
        linProbs = new double[0];// initialize the tree and rates interval counter
        int treeInterval = 0, ratesInterval = 0;        
        double nextEventTime = 0.0;
		coalescentRates = mascotInput.get().dynamicsInput.get().getCoalescentRate(ratesInterval);  
        migrationRates = mascotInput.get().dynamicsInput.get().getBackwardsMigration(ratesInterval);
		indicators = mascotInput.get().dynamicsInput.get().getIndicators(ratesInterval);  
		        
        // Time to the next rate shift or event on the tree
        double nextTreeEvent = mascotInput.get().treeIntervalsInput.get().getInterval(treeInterval);
        double nextRateShift = mascotInput.get().dynamicsInput.get().getInterval(ratesInterval);
        double collapsTime = mascotInput.get().collapsTimeInput.get().getValue();
        
        double time = 0.0;
        // Calculate the likelihood
        do {       
        	nextEventTime = Math.min(Math.min(nextTreeEvent, nextRateShift), collapsTime); 	 
        	time += nextEventTime;
        	if (nextEventTime > 0) {													// if true, calculate the interval contribution        		
                if(recalculateLogP){
                	System.out.println(logP);
    				System.err.println("ode calculation stuck, reducing tolerance, new tolerance= " + maxTolerance);
    				maxTolerance *=0.9;
    		    	recalculateLogP = false;
    				System.exit(0);
                	return calculateNodeStates();
                }
                if(stepSizeInput.get()!=null){
		        	double[] linProbs_for_ode = new double[linProbs.length]; 
	                FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(stepSizeInput.get());	                
//	                integrator.setMaxEvaluations((int) 1e5);  // set the maximal number of evaluations              
	                FirstOrderDifferentialEquations ode = new MascotODE(migrationRates, coalescentRates, nrLineages , coalescentRates.length);                
	                try {
	                	integrator.integrate(ode, 0, linProbs, nextEventTime, linProbs_for_ode);
	                }catch(Exception e){
	                	System.out.println(e);
	                	recalculateLogP = true;    				
	                }                
		            for (int i = 0; i<linProbs.length; i++)
		        		linProbs[i] = linProbs_for_ode[i]; 
	        	}else {
	        		if (coalescentRates.length>1){
		        		Euler2ndOrder euler;
		        		if (mascotInput.get().dynamicsInput.get().hasIndicators)
		        			euler = new Euler2ndOrder(migrationRates, indicators, coalescentRates, nrLineages , coalescentRates.length, epsilonInput.get(), maxStepInput.get());
		        		else
		        			euler = new Euler2ndOrder(migrationRates, coalescentRates, nrLineages , coalescentRates.length, epsilonInput.get(), maxStepInput.get());
		        		
			        	double[] linProbs_tmp = new double[linProbs.length+1]; 
			        	double[] linProbs_tmpdt = new double[linProbs.length+1]; 
			        	double[] linProbs_tmpddt = new double[linProbs.length+1]; 
			        	double[] linProbs_tmpdddt = new double[linProbs.length+1]; 
			        	
			        	for (int i = 0; i < linProbs.length; i++) linProbs_tmp[i] = linProbs[i];		        	
			        	
			        	linProbs[linProbs.length-1] = 0;
			        	euler.calculateValues(nextEventTime, linProbs_tmp, linProbs_tmpdt, linProbs_tmpddt, linProbs_tmpdddt, linProbs.length + 1);		        	
		        		
			            for (int i = 0; i < linProbs.length; i++) linProbs[i] = linProbs_tmp[i]; 
			            
			            logP += linProbs_tmp[linProbs_tmp.length-1];
//	        			System.out.println(time);
	        		}else{
	        			logP -= nextEventTime*coalescentRates[0]*activeLineages.size()*(activeLineages.size()-1)/2;
	        		}
	        	}        	
        	}
       	
        	if (nextTreeEvent <= Math.min(nextRateShift, collapsTime)){
 	        	if (mascotInput.get().treeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.COALESCENT) {
// 	        		System.out.print(String.format("%.3f ", nextTreeEvent));
	        		logP += normalizeLineages();									// normalize all lineages before event		
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		logP += coalesce(treeInterval, ratesInterval);	  				// calculate the likelihood of the coalescent event
	        	}
 	       		
 	       		if (mascotInput.get().treeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	       			
 	       			if (linProbs.length > 0)
 	       				logP += normalizeLineages();								// normalize all lineages before event
	       			nrLineages++;													// sampling event increases the number of lineages by one
 	       			sample(treeInterval, ratesInterval);							// calculate the likelihood of the sampling event if sampling rate is given
	       		}	
 	       		
 	       		treeInterval++;
        		nextRateShift -= nextTreeEvent;   
        		collapsTime -= nextTreeEvent;  
        		try{
        			nextTreeEvent = mascotInput.get().treeIntervalsInput.get().getInterval(treeInterval);
        		}catch(Exception e){
        			break;
        		}
        	}else if(nextRateShift < collapsTime){
        		System.err.println("rate shift are not allowed");
        		System.exit(0);
        	}else{
        		coalescentRates = new double[1];
        		migrationRates = new double[1];
        		coalescentRates[0] = 1/(2*mascotInput.get().ancestralNeInput.get().getValue());
                migrationRates[0] = 0.0;
        		nextTreeEvent -= collapsTime;
        		nextRateShift -= collapsTime;
 	       		collapsTime = Double.POSITIVE_INFINITY;
        	}
        	if (logP == Double.NEGATIVE_INFINITY)
        		return logP;
        }while(nextTreeEvent <= Double.POSITIVE_INFINITY);
//        System.exit(0);
//        System.out.print("\n");
//        first = false;
//        System.exit(0);
        
//        System.out.println(times[1]/(times[0]+times[1]+times[2]+times[3]+times[4]));
//        System.out.println(times[0]/(times[0]+times[1]) + " " + times[2]/(times[2]+times[3]));
        
//        System.out.println((times[0]/(times[0]+times[1]+times[2]+times[3])) + " " + times[1]/(times[0]+times[1]+times[2]+times[3]) +
//        		" " + times[2]/(times[0]+times[1]+times[2]+times[3]) + " " + times[3]/(times[0]+times[1]+times[2]+times[3]));
		return logP;  	
    }   
    
    
    private double normalizeLineages(){
    	if (linProbs==null)
    		return 0.0;
    		
    	if (coalescentRates.length==1)
    		return 0.0;
    	
    	double interval = 0.0;
    	for (int i = 0; i < nrLineages; i++){
    		double lineProbs = 0.0;
    		for (int j = 0; j < states; j++)
    			if (linProbs[i*states+j]>=0.0){
    				lineProbs += linProbs[i*states+j];
    			}else{
    				// try recalculation after lowering the tolerance
    				System.out.println("dkfslj");
    				System.out.println(linProbs[i*states+j]);
    				recalculateLogP = true;
    				return Math.log(1.0);
    			}
    		for (int j = 0; j < states; j++){
    			if (lineProbs==0.0)
    				return Double.NEGATIVE_INFINITY;
    			linProbs[i*states+j] = linProbs[i*states+j]/lineProbs;
    		}    		
    		interval +=lineProbs;
    	}	
		// return mean P_t(T)
		return Math.log(interval/(nrLineages));

    }
    
    private void sample(int currTreeInterval, int currRatesInterval) {
		List<Integer> incomingLines = mascotInput.get().treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		int newLength = linProbs.length + incomingLines.size()*states;
		
		double[] linProbsNew = new double[newLength];
		
		
		for (int i = 0; i < linProbs.length; i++)
			linProbsNew[i] = linProbs[i];
		
		
		int currPosition = linProbs.length;
		
		/*
		 * If there is no trait given as Input, the model will simply assume that
		 * the last value of the taxon name, the last value after a _, is an integer
		 * that gives the type of that taxon
		 */
		if (mascotInput.get().dynamicsInput.get().typeTraitInput.get()!=null){
			for (Integer l : incomingLines) {
				activeLineages.add(l);//.getNr());
				int sampleState = mascotInput.get().dynamicsInput.get().getValue(tree.getNode(l).getID());
				
				if (sampleState>= mascotInput.get().dynamicsInput.get().getDimension()){
					System.err.println("sample discovered with higher state than dimension");
//					System.exit(1);
				}
				
				for (int i = 0; i < states; i++){
					if (i == sampleState){
						linProbsNew[currPosition] = 1.0;currPosition++;
					}
					else{
						linProbsNew[currPosition] = 0.0;currPosition++;
					}
				}
			}				
		}else{
			for (Integer l : incomingLines) {
				activeLineages.add(l);//.getNr());
				String sampleID = tree.getNode(l).getID();
				int sampleState = 0;
				if (states > 1){				
					String[] splits = sampleID.split("_");
					sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
				}
				for (int i = 0; i < states; i++){
					if (i == sampleState){
						linProbsNew[currPosition] = 1.0;currPosition++;
					}
					else{
						linProbsNew[currPosition] = 0.0;currPosition++;
					}
				}
			}	
		}
		linProbs = linProbsNew;
		// store the node
//		storeNode(currTreeInterval, currRatesInterval, linProbs, logP, activeLineages);
    }
    
    private double coalesce(int currTreeInterval, int currRatesInterval) {
		List<Integer> coalLines = mascotInput.get().treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
    	if (coalLines.size() > 2) {
			System.err.println("Unsupported coalescent at non-binary node");
			System.exit(0);
		}
    	if (coalLines.size() < 2) {
    		System.out.println();
    		System.out.println("WARNING: Less than two lineages found at coalescent event!");
    		System.out.println();
    		return Double.NaN;
		}
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0));//.getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1));//.getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println(coalLines.get(0)/*.getNr()*/ + " " + coalLines.get(1)/*.getNr()*/ + " " + activeLineages);
			System.out.println("daughter lineages at coalescent event not found");
			return Double.NaN;
		}
		DoubleMatrix lambda = DoubleMatrix.zeros(states);
		
		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */
		if (coalescentRates.length>1){
	        for (int k = 0; k < states; k++) { 
	        	Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];			
				if (!Double.isNaN(pairCoalRate)){
					lambda.put(k, pairCoalRate);
				}else{
					return Double.NEGATIVE_INFINITY;
				}
	        }
		}else{
			for (int i = 0; i < lambda.length; i++)
				lambda.put(i, 1/states);			
		}
        
        activeLineages.add(tree.getNode(coalLines.get(0)).getParent().getNr());        
        
        // get the node state probabilities
		DoubleMatrix pVec = new DoubleMatrix();
		pVec.copy(lambda);
		pVec = pVec.div(pVec.sum());
				
		stateProbabilities[tree.getNode(coalLines.get(0)).getParent().getNr() - nrSamples] = pVec;
		
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
			linProbsNew[linCount*states + j] = pVec.get(j);
		}
		// set p to pnew
		linProbs = linProbsNew;	
		
		
		//Remove daughter lineages from the line state probs
		if (daughterIndex1>daughterIndex2){
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex1);
			activeLineages.remove(daughterIndex2);			
		}else{
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex2);
			activeLineages.remove(daughterIndex1);			
		}		
     
		if (lambda.min()<0.0){
			System.err.println("Coalescent probability is: " + lambda.min());
			return Double.NEGATIVE_INFINITY;
		}				
				
		leftID[tree.getNode(coalLines.get(0)).getParent().getNr() - nrSamples] = tree.getNode(coalLines.get(0)).getParent().getLeft().getNr();
		rightID[tree.getNode(coalLines.get(0)).getParent().getNr() - nrSamples] = tree.getNode(coalLines.get(0)).getParent().getRight().getNr();	
		
		if (coalescentRates.length>1){
			if (lambda.sum()==0)
				return Double.NEGATIVE_INFINITY;
			else
				return Math.log(lambda.sum());
		}else{
			return Math.log(coalescentRates[0]);						
		}
    }
            

    
	private DoubleMatrix getStateProb(int nr) {
		used[nr - nrSamples] = true;
		return stateProbabilities[nr - nrSamples];
	}


	
}
