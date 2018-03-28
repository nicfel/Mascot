package beast.mascot.logger;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
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
import beast.evolution.tree.coalescent.IntervalType;
import beast.mascot.distribution.Mascot;
import beast.mascot.ode.Euler2ndOrderTransitions;
import beast.mascot.ode.MascotODEUpDown;


/**
 * @author Nicola Felix Mueller (nicola.felix.mueller@gmail.com)
 */
public class StructuredTreeLogger extends Tree implements Loggable {
	
	
	public Input<Mascot> mascotInput = new Input<>("mascot", "Input of rates", Input.Validate.REQUIRED);
    
	public Input<Double> epsilonInput = new Input<>("epsilon", "step size for the RK4 integration",0.00001);
	public Input<Double> maxStepInput = new Input<>("maxStep", "step size for the RK4 integration", Double.POSITIVE_INFINITY);
	public Input<Double> stepSizeInput = new Input<>("stepSize", "step size for the RK4 integration");
	
	public Input<Boolean> useUpDown = new Input<>("upDown", "if up down algorithm is to use for the node state calculation", true);
    
    
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
    
    boolean updown = true;

    private DecimalFormat df;
    private String type;
    
    private int states;
    boolean[] used;
	boolean report;
    

	
    @Override
    public void initAndValidate() {  	
    	
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
		        
		        for (int i = 0 ; i < states; i++)
		        	buf.append(String.format("%s=%.3f,", mascotInput.get().dynamicsInput.get().getStringStateValue(i), stateProbs.get(i)));
		        
//		        buf.append(String.format("%.3f", stateProbs.get(states-1)));
//		        buf.append("}");
	        
		        buf.append("max=");
		        buf.append(String.format("%s", 
		        		mascotInput.get().dynamicsInput.get().getStringStateValue(stateProbs.argmax())));
		        
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

	
    private void CalculateNodeStates() throws Exception{  
    	// newly calculate tree intervals
    	mascotInput.get().treeIntervalsInput.get().calculateIntervals();
    	// correctly calculate the daughter nodes at coalescent intervals in the case of
    	// bifurcation or in case two nodes are at the same height
    	mascotInput.get().treeIntervalsInput.get().swap();    	
 
    	leftID = new int[mascotInput.get().treeIntervalsInput.get().getSampleCount()];
    	rightID = new int[mascotInput.get().treeIntervalsInput.get().getSampleCount()];

    	
    	stateProbabilities = new DoubleMatrix[mascotInput.get().treeIntervalsInput.get().getSampleCount()];
    	stateProbabilitiesDown = new DoubleMatrix[mascotInput.get().treeIntervalsInput.get().getSampleCount()];
    	TransitionProbabilities = new DoubleMatrix[mascotInput.get().treeIntervalsInput.get().getSampleCount()*2];        
        nrSamples = mascotInput.get().treeIntervalsInput.get().getSampleCount() + 1;        

    	
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages = new ArrayList<Integer>(); 

        nrLineages = 0;
        linProbs = new double[0];// initialize the tree and rates interval counter
        transitionProbs = new double[0];// initialize the tree and rates interval counter
        
        int treeInterval = 0, ratesInterval = 0;        
        double nextEventTime = 0.0;
		coalescentRates = mascotInput.get().dynamicsInput.get().getCoalescentRate(ratesInterval);  
        migrationRates = mascotInput.get().dynamicsInput.get().getBackwardsMigration(ratesInterval);
		indicators = mascotInput.get().dynamicsInput.get().getIndicators(ratesInterval);  
        // Time to the next rate shift or event on the tree
        double nextTreeEvent = mascotInput.get().treeIntervalsInput.get().getInterval(treeInterval);
        double nextRateShift = mascotInput.get().dynamicsInput.get().getInterval(ratesInterval);

             

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
	        		if (mascotInput.get().dynamicsInput.get().hasIndicators)
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
 	        	if (mascotInput.get().treeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.COALESCENT) {
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		normalizeLineages();									// normalize all lineages before event		
        			coalesce(treeInterval);	  				// calculate the likelihood of the coalescent event
	        	}
 	       		
 	       		if (mascotInput.get().treeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	
 	       			nrLineages++;													// sampling event increases the number of lineages by one
 	       			if (linProbs.length > 0)
 	       				normalizeLineages();								// normalize all lineages before event
 	       			sample(treeInterval);							// calculate the likelihood of the sampling event if sampling rate is given
	       		}	
 	       		
 	       		treeInterval++;
        		nextRateShift -= nextTreeEvent;   
        		try{
        			nextTreeEvent = mascotInput.get().treeIntervalsInput.get().getInterval(treeInterval);
        		}catch(Exception e){
        			break;
        		}
        	}else{
        		ratesInterval++;
        		coalescentRates = mascotInput.get().dynamicsInput.get().getCoalescentRate(ratesInterval);  
                migrationRates = mascotInput.get().dynamicsInput.get().getBackwardsMigration(ratesInterval);
        		indicators = mascotInput.get().dynamicsInput.get().getIndicators(ratesInterval);  
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = mascotInput.get().dynamicsInput.get().getInterval(ratesInterval);
        	}
        	
        }while(nextTreeEvent <= Double.POSITIVE_INFINITY);
        currTreeInterval = mascotInput.get().treeIntervalsInput.get().getIntervalCount()-1;
        
        do{
		  	if (mascotInput.get().treeIntervalsInput.get().getIntervalType(currTreeInterval) == IntervalType.COALESCENT) {
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
		List<Node> incomingLines = mascotInput.get().treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		// calculate the new length of the arrays for the transition and lineage states
		int newLengthLineages = linProbs.length + incomingLines.size()*states;
		int newLengthTransitions = transitionProbs.length + incomingLines.size()*states*states;
		
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
		if (mascotInput.get().dynamicsInput.get().typeTraitInput.get()!=null){
			for (Node l : incomingLines) {
				activeLineages.add(l.getNr());
				int sampleState = (int) mascotInput.get().dynamicsInput.get().getValue(l.getID());
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
			for (Node l : incomingLines) {
				activeLineages.add(l.getNr());
				String sampleID = l.getID();
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

    	

		List<Node> coalLines = mascotInput.get().treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
    	if (coalLines.size() > 2) {
			System.err.println("Unsupported coalescent at non-binary node");
			System.exit(0);
		}
    	if (coalLines.size() < 2) {
    		System.out.println();
    		System.out.println("WARNING: Less than two lineages found at coalescent event!");
    		System.exit(0);
		}
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println("daughter lineages at coalescent event not found");
    		System.exit(0);
		}
		DoubleMatrix lambda = DoubleMatrix.zeros(states);

		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */
        for (int k = 0; k < states; k++) { 
        	Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];			
			if (!Double.isNaN(pairCoalRate)){
				lambda.put(k, pairCoalRate);
			}else{
//	    		System.exit(0);
			}
        }
        
        activeLineages.add(coalLines.get(0).getParent().getNr());
        
        
        // get the node state probabilities
		DoubleMatrix pVec = new DoubleMatrix();
		pVec.copy(lambda);
		pVec = pVec.div(pVec.sum());
		
		// save the node states conditioned on the subtree
		stateProbabilities[coalLines.get(0).getParent().getNr() - nrSamples] = pVec;

		
		// get the transition probabilities of daughter lineage 1
		DoubleMatrix tP1 = DoubleMatrix.zeros(states,states);
		for (int i = 0; i< states; i++){
			for (int j = 0; j< states; j++){
				tP1.put(i, j, transitionProbs[daughterIndex1*states*states+i*states+j]);
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
		DoubleMatrix tP2 = DoubleMatrix.zeros(states,states);
		for (int i = 0; i< states; i++){
			for (int j = 0; j< states; j++){
				tP2.put(i, j, transitionProbs[daughterIndex2*states*states+i*states+j]);
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
			linProbsNew[linCount*states + j] = pVec.get(j);
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
		TransitionProbabilities[coalLines.get(0).getNr()] = tP1;
		TransitionProbabilities[coalLines.get(1).getNr()] = tP2;
		
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
		
		if(coalLines.get(0).getParent().getNr() != mascotInput.get().treeIntervalsInput.get().treeInput.get().getNode(coalLines.get(1).getNr()).getParent().getNr())
			System.err.println("wrong daughter parent");
		if(coalLines.get(1).getParent().getNr() != mascotInput.get().treeIntervalsInput.get().treeInput.get().getNode(coalLines.get(0).getNr()).getParent().getNr())
			System.err.println("wrong daughter parent");
		if(coalLines.get(1).getParent().getNr() != coalLines.get(0).getParent().getNr())
			System.err.println("coalescent nodes don't have the same parent");
		
		leftID[coalLines.get(0).getParent().getNr() - nrSamples] = coalLines.get(0).getParent().getLeft().getNr();
		rightID[coalLines.get(0).getParent().getNr() - nrSamples] = coalLines.get(0).getParent().getRight().getNr();
    }

    private void coalesceDown(int currTreeInterval) {
		List<Node> parentLines = mascotInput.get().treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		if (parentLines.size()!=1){
			System.err.println("to many lineages, while coalescening down");
			System.exit(0);
		}
		Node parentNode = parentLines.get(0);
		
		if (!parentNode.isRoot()){
			DoubleMatrix start = stateProbabilities[parentNode.getNr() - nrSamples];
			DoubleMatrix end = stateProbabilitiesDown[parentNode.getParent().getNr() - nrSamples];
			DoubleMatrix flow = TransitionProbabilities[parentNode.getNr()];
			DoubleMatrix otherSideInfo = end.div(start.transpose().mmul(flow));
			// get rid of NaN from division by 0
			for (int i = 0; i < otherSideInfo.length; i++)
				if (Double.isNaN(otherSideInfo.get(i)))
					otherSideInfo.put(i, 0.0);
				
				
				
			DoubleMatrix conditional = flow.mmul(otherSideInfo);
			conditional = conditional.mul(start);
			stateProbabilitiesDown[parentNode.getNr() - nrSamples] = conditional.div(conditional.sum());
//			if (!(conditional.get(0) >= 0.0 && conditional.get(0)<=1.0))
//				conditional.print();
		}else{
//			DoubleMatrix d1 = stateProbabilities[parentNode.getLeft().getNr() - nrSamples];
//			DoubleMatrix d2 = stateProbabilities[parentNode.getRight().getNr() - nrSamples];
//			DoubleMatrix f1 = TransitionProbabilities[parentNode.getLeft().getNr()];
//			DoubleMatrix f2 = TransitionProbabilities[parentNode.getLeft().getNr()];
//			d1.transpose().mmul(f1).print();
//			d2.transpose().mmul(f2).print();
//			d1.print();
//			d2.print();
//			System.out.println();
//			f1.print();
//			f2.print();
//			System.out.println();
//			stateProbabilities[parentNode.getNr() - nrSamples].print();
//			
//			System.out.println();
//			System.exit(0);
			stateProbabilitiesDown[parentNode.getNr() - nrSamples] = stateProbabilities[parentNode.getNr() - nrSamples];
    	}
	}

    
	private DoubleMatrix getStateProb(int nr) {
		if(useUpDown.get()){			
			used[nr - nrSamples] = true;
			return stateProbabilitiesDown[nr - nrSamples] ;
		}else{
			used[nr - nrSamples] = true;
			return stateProbabilities[nr - nrSamples] ;
		}
	}


	
}
