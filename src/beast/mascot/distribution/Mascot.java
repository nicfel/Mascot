package beast.mascot.distribution;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.jblas.DoubleMatrix;

import beast.core.CalculationNode;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.coalescent.IntervalType;
import beast.mascot.dynamics.Dynamics;
import beast.mascot.ode.*;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
@Citation("Nicola F. Müller, David A. Rasmussen, Tanja Stadler (2017)\n  The Structured Coalescent and its Approximations.\n  Mol Biol Evol 2017 msx186. doi: 10.1093/molbev/msx186")
public class Mascot extends StructuredTreeDistribution {
	
	public Input<Dynamics> dynamicsInput = new Input<>("dynamics", "Input of rates", Input.Validate.REQUIRED);
	public Input<Double> epsilonInput = new Input<>("epsilon", "step size for the RK4 integration",0.001);
	public Input<Double> maxStepInput = new Input<>("maxStep", "step size for the RK4 integration", Double.POSITIVE_INFINITY);
//	public Input<Double> stepSizeInput = new Input<>("stepSize", "step size for the RK4 integration");
	public Input<Boolean> saveRamInput = new Input<>("saveRamInput", "doesn't save intermediate steps", false);
	
	final static int N = 10; // cache size
    
	public int samples;
	public int nrSamples;
	public DoubleMatrix[] stateProbabilities;
    
    private int nrLineages;   

    // current rates         
    private double[] migrationRates;
    private int [] indicators;
    private double[] coalescentRates; 	

    
    // Set up for lineage state probabilities
    private ArrayList<Integer> activeLineages;
	private double[] linProbs;
	private double[] linProbsNew;
	private int linProbsLength;
	private int states;
	
    // store the linProbs, multiplicators and logP's at coalescent points in jagged arrays from last time
    private double[][] coalLinProbs;
    private double[] coalLogP;
    private int[] coalRatesInterval;
    private ArrayList<Integer> [] coalActiveLineages;
	private double [] nextTreeEvents;
	private double [] nextRateShifts;
    
    // deep store the things above for MCMC
    private double[][] storeLinProbs;
    private double[] storeLogP;
    private int[] storeRatesInterval;
    private ArrayList<Integer> [] storeActiveLineages;
	private double [] storedNextTreeEvents;
	private double [] storedNextRateShifts;

	private boolean storedFirst;
	//private int [] treeIntervalNrs;
	//private int [] storedTreeIntervalNrs;
	//private int [] lineagesAddded;
	//private int [] storedLineagesAddded;

    // check if this is the first calculation
    private int first = 0;

	
	// maximum integration error tolerance
    private double maxTolerance = 1e-3;            
    private boolean recalculateLogP;
	Euler2ndOrderBase euler;
	TreeInterface tree;
	Dynamics dynamics;
	StructuredTreeIntervals treeIntervals;

	int [] nodeType;

    @Override
    public void initAndValidate(){
    	dynamics = dynamicsInput.get();
    	treeIntervals = treeIntervalsInput.get();
    	tree = treeInput.get();
    	if (tree == null) {
    		tree = treeIntervals.treeInput.get();
    	}
    	treeIntervals.calculateIntervals();       
    	stateProbabilities = new DoubleMatrix[treeIntervals.getSampleCount()];
        nrSamples = treeIntervals.getSampleCount() + 1;    
        states = dynamics.getDimension();
                

    	// initialize storing arrays and ArrayLists
    	coalLinProbs = new double[N][];
    	storeLinProbs = new double[N][];

    	coalLogP = new double[N];
    	storeLogP = new double[N];
    	
    	coalRatesInterval = new int[N];
    	storeRatesInterval = new int[N];
    	
    	coalActiveLineages = new ArrayList[N];
    	storeActiveLineages = new ArrayList[N];
    	
    	nextTreeEvents = new double[N];
    	storedNextTreeEvents = new double[N];

    	nextRateShifts = new double[N];
    	storedNextRateShifts = new double[N];

    	for (int i = 0; i < N; i++) {
    		coalActiveLineages[i] = new ArrayList<>();;
    		storeActiveLineages[i] = new ArrayList<>();;
    	}

    	activeLineages = new ArrayList<>();
    	
    	int intCount = treeIntervals.getIntervalCount();
    	int MAX_SIZE = intCount * states;
    	linProbs_for_ode = new double[MAX_SIZE];
    	linProbs_tmp = new double[MAX_SIZE];
    	linProbs = new double[MAX_SIZE];
    	linProbsNew = new double[MAX_SIZE];
    	
    	if (Euler2ndOrderNative.loadLibrary()) {
    		euler = new Euler2ndOrderNative();
    	} else {
    		switch (states) {
    		case 2: euler = new Euler2ndOrder2(); break;
    		case 3: euler = new Euler2ndOrder3(); break;
    		case 4: euler = new Euler2ndOrder4(); break;
    		case 5: euler = new Euler2ndOrder5(); break;
    		case 6: euler = new Euler2ndOrder6(); break;
    		case 7: euler = new Euler2ndOrder7(); break;
    		case 8: euler = new Euler2ndOrder8(); break;
    		case 9: euler = new Euler2ndOrder9(); break;
    		case 10: euler = new Euler2ndOrder10(); break;
    		default: euler = new Euler2ndOrder(); break;
    		}
    	}
    	Log.warning("Using " + euler.getClass().getSimpleName());
    	
    	euler.setup(MAX_SIZE, states, epsilonInput.get(), maxStepInput.get());
    	
    	nodeType = new int[tree.getNodeCount()];
    	for (int i = 0; i < tree.getLeafNodeCount(); i++) {
    		nodeType[i] = dynamics.getValue(tree.getNode(i).getID());
    	}
    }
    
    double [] linProbs_for_ode;
    double [] linProbs_tmp;

    public double calculateLogP() {    	
    	// newly calculate tree intervals
    	treeIntervals.calculateIntervals();
    	// correctly calculate the daughter nodes at coalescent intervals in the case of
    	// bifurcation or in case two nodes are at the same height
    	treeIntervals.swap();    	
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages.clear();
        logP = 0;
        nrLineages = 0;
        //linProbs = new double[0];// initialize the tree and rates interval counter
        linProbsLength = 0;
        int treeInterval = 0, ratesInterval = 0;        
     	double nextEventTime = 0.0;
        
        // Time to the next rate shift or event on the tree
        double nextTreeEvent = treeIntervals.getInterval(treeInterval);
        double nextRateShift = dynamics.getInterval(ratesInterval);
        
      //ArrayList<Integer> oldActiveLineages = new ArrayList<>();
        
      if (first > 0 && !dynamics.isDirtyCalculation() 
				&& treeIntervals.firstDirtyInterval > 2) {
    	  
			treeInterval = treeIntervals.firstDirtyInterval - 1;
			ratesInterval = restoreNode(treeInterval-1);
			nextTreeEvent = nextTreeEvents[treeInterval - 1];
			nextRateShift = nextRateShifts[treeInterval - 1];
  
//    	  
//    	  
//        // restore the likelihood to last known good place
//        do {
//        	//oldActiveLineages.clear();
//        	//oldActiveLineages.addAll(activeLineages);
//        	
//        	nextEventTime = Math.min(nextTreeEvent, nextRateShift); 	 
//
//    		if (nextTreeEvent <= nextRateShift){
//        		// Check if the last interval was reached
//        		if (treeInterval == treeIntervals.intervalCount){
//        			Log.warning("Restoring to the finish!");
//        			logP = coalLogP[coalLogP.length-1];
//        			return logP;
//        		}
//        		boolean isDirty = false;
//        		if (treeIntervals.getIntervalType(treeInterval) == IntervalType.COALESCENT) {        		
////	        		logP += normalizeLineages(linProbs);							// normalize all lineages before event		
// 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
// 	        		isDirty = coalesceX(treeInterval, ratesInterval);	  				// calculate the likelihood of the coalescent event
//	        	}
// 	       		
// 	       		if (treeIntervals.getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	       			
//// 	       			if (linProbsLength > 0)
//// 	       				logP += normalizeLineages(linProbs);						// normalize all lineages before event
//	       			nrLineages++;													// sampling event increases the number of lineages by one
//	       			isDirty = sampleX(treeInterval, ratesInterval);							// calculate the likelihood of the sampling event if sampling rate is given
//	       		}	
// 	       		
// 	       	
// 	    		if (treeIntervals.intervalIsDirty(treeInterval+1) || isDirty) { //differ(activeLineages, coalActiveLineages.get(treeInterval))) {
// 	    			if (treeInterval <= 2) {
// 	        			//Log.warning("Reset @ " + first);
// 	        			treeInterval = 0;
// 	        			ratesInterval = 0;
// 	        	        activeLineages = new ArrayList<Integer>();
// 	        	        logP = 0;
// 	        	        nrLineages = 0;
// 	        	        linProbsLength = 0;
// 	        	     	nextEventTime = 0.0;
// 	        	        nextTreeEvent = treeIntervals.getInterval(treeInterval);
// 	        	        nextRateShift = dynamics.getInterval(ratesInterval);        			
// 	        			break; 	    				
// 	    			}
// 	    			
// 	    			
//// 	    			Log.warning.println("Restore I " + treeInterval + " " + treeIntervals.firstDirtyInterval);
// 	    			treeInterval++;
// 	    			ratesInterval = restoreNode(treeInterval-2);
// 	    			//activeLineages = oldActiveLineages;
//        			//nextTreeEvent = treeIntervals.getInterval(treeInterval);
// 	    			nextTreeEvent = nextTreeEvents[treeInterval-1];
// 	    			nextRateShift = nextRateShifts[treeInterval-1];
// 	    			treeInterval--;
// 	        		//nextRateShift -= nextTreeEvent;   
// 	 	       		//treeInterval++;
// 	    			break;
// 	    		}
// 	       		treeInterval++;
//        		nextRateShift -= nextTreeEvent;   
//        		try{
//        			nextTreeEvent = treeIntervals.getInterval(treeInterval);
//        		} catch(Exception e) {
//        			break;
//        		}
//        	} else {
//        		ratesInterval++;
//        		coalescentRates = dynamics.getCoalescentRate(ratesInterval);  
//                migrationRates = dynamics.getBackwardsMigration(ratesInterval);
//        		indicators = dynamics.getIndicators(ratesInterval);  
//        		nextTreeEvent -= nextRateShift;
// 	       		nextRateShift = dynamics.getInterval(ratesInterval);
//        	}
//        	//if (logP == Double.NEGATIVE_INFINITY) {
//        	//	return logP;
//        	//}
//        } while(nextTreeEvent <= Double.POSITIVE_INFINITY);
    }


		coalescentRates = dynamics.getCoalescentRate(ratesInterval);  
        migrationRates = dynamics.getBackwardsMigration(ratesInterval);
		indicators = dynamics.getIndicators(ratesInterval);
		nrLineages = activeLineages.size();

        // Calculate the likelihood
        do {       
        	nextEventTime = Math.min(nextTreeEvent, nextRateShift); 	 
        	if (treeInterval >= 235 && first == 5) {
        		int h = 3;
        		h++;
        	}
       	
        	if (nextEventTime > 0) {													// if true, calculate the interval contribution        		
                if (recalculateLogP) {
    				System.err.println("ode calculation stuck, reducing tolerance, new tolerance= " + maxTolerance);
    				maxTolerance *=0.9;
    		    	recalculateLogP = false;
    				System.exit(0);
                	return calculateLogP();
                }
        		logP += doEuler(nextEventTime);
        	}
       	
        	if (nextTreeEvent <= nextRateShift){
 	        	if (treeIntervals.getIntervalType(treeInterval) == IntervalType.COALESCENT) {
// 	        		System.out.print(String.format("%.3f ", nextTreeEvent));
	        		logP += normalizeLineages(linProbs);									// normalize all lineages before event		
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		logP += coalesce(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);	  				// calculate the likelihood of the coalescent event
	        	}
 	       		
 	       		if (treeIntervals.getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	       			
 	       			//if (linProbsLength > 0)
 	       			//	logP += normalizeLineages(linProbs);	     obsolete: doEuler already takes care of normalisation
	       			nrLineages++;													// sampling event increases the number of lineages by one
 	       			sample(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);							// calculate the likelihood of the sampling event if sampling rate is given
	       		}	
 	       		
 	       		treeInterval++;
        		nextRateShift -= nextTreeEvent;   
        		try{
        			nextTreeEvent = treeIntervals.getInterval(treeInterval);
        		}catch(Exception e){
        			break;
        		}
        	} else {
        		ratesInterval++;
        		coalescentRates = dynamics.getCoalescentRate(ratesInterval);  
                migrationRates = dynamics.getBackwardsMigration(ratesInterval);
        		indicators = dynamics.getIndicators(ratesInterval);  
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = dynamics.getInterval(ratesInterval);
        	}
        	if (logP == Double.NEGATIVE_INFINITY) {
        		return logP;
        	}
        } while(nextTreeEvent <= Double.POSITIVE_INFINITY);

        first++;
//        System.err.println(logP);
		return logP;  	
    }

    private boolean differ(ArrayList<Integer> a1, ArrayList<Integer> a2) {
		if (a1.size() != a2.size()) {
			return true;
		}
		for (int i = 0; i < a2.size(); i++) {
			if (a1.get(i) - a2.get(i) != 0) {
				return true;
			}
		}
		return false;
	}

	double [] storedMigrationRates = new double[0];
    double [] storedCoalescentRates = new double[0];
    int storedNrLineages = -1;
    
	private double doEuler(double nextEventTime) {
		//for (int i = 0; i < linProbs.length; i++) linProbs_tmp[i] = linProbs[i];
		if (linProbs_tmp.length != linProbsLength + 1) {
			linProbs_tmp= new double[linProbsLength + 1];
		}
		System.arraycopy(linProbs,0,linProbs_tmp,0,linProbsLength);
		linProbs_tmp[linProbsLength] = 0;

		linProbs[linProbsLength-1] = 0;
		

		if (dynamics.hasIndicators) {
			euler.initWithIndicators(migrationRates, indicators, coalescentRates, nrLineages);
			euler.calculateValues(nextEventTime, linProbs_tmp, linProbsLength + 1);
		} else {
			euler.initAndcalculateValues(migrationRates, coalescentRates, nrLineages, nextEventTime, linProbs_tmp, linProbsLength + 1);
		}
		
		//		System.out.println(Arrays.toString(linProbs));		

		//for (int i = 0; i < linProbs.length; i++) linProbs[i] = linProbs_tmp[i];
		System.arraycopy(linProbs_tmp,0,linProbs,0,linProbsLength);

		return linProbs_tmp[linProbsLength];
	}

	
//    private void integrate(double duration){
//    }
//
//    private void ei(double duration, double[] linProbs_for_ode, double[] meanLinProbs){
//    	eulerIntegration(duration, linProbs_for_ode, meanLinProbs);   	
//    }
    
    private double normalizeLineages(double [] linProbs){
    	if (linProbs==null)
    		return 0.0;
    	
    	
    	double interval = 0.0;
    	for (int i = 0; i < nrLineages; i++) {
    		double lineProbs = 0.0;
    		int u = i * states;
    		for (int j = 0; j < states; j++) {
    			if (linProbs[u]>=0.0){
    				lineProbs += linProbs[u];
    			} else {
    				// try recalculation after lowering the tolerance
    				System.out.println(linProbs[u]);
    				recalculateLogP = true;
    				return Math.log(1.0);
    			}
    			u++;
    		}
			if (lineProbs==0.0) {
				return Double.NEGATIVE_INFINITY;
			}
			u = i * states;
    		for (int j = 0; j < states; j++) {
    			linProbs[u] = linProbs[u]/lineProbs;
    			u++;
    		}    		
    		interval += lineProbs;
    	}	
		// return mean P_t(T)
		return Math.log(interval/(nrLineages));

    }
    
    private void sample(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
		int incomingLines = treeIntervals.getLineagesAdded(currTreeInterval);
		int newLength = linProbsLength + 1 * states;
		
		//double [] linProbsNew = new double[newLength];
//		System.arraycopy(linProbs, 0, linProbsNew, 0, linProbsLength);
//		for (int i = 0; i < linProbs.length; i++)
//			linProbsNew[i] = linProbs[i];
		
		
		int currPosition = linProbsLength;
		
		/*
		 * If there is no trait given as Input, the model will simply assume that
		 * the last value of the taxon name, the last value after a _, is an integer
		 * that gives the type of that taxon
		 */
		if (dynamics.typeTraitInput.get()!=null){
			Integer l = incomingLines; {
				activeLineages.add(l);//.getNr());
				int sampleState = nodeType[l];//dynamics.getValue(tree.getNode(l).getID());
				
				if (sampleState>= dynamics.getDimension()){
					System.err.println("sample discovered with higher state than dimension");
//					System.exit(1);
				}
				
				for (int i = 0; i < states; i++){
					if (i == sampleState){
						linProbs[currPosition] = 1.0;currPosition++;
					}
					else{
						linProbs[currPosition] = 0.0;currPosition++;
					}
				}
			}				
		}else{
			Integer l = incomingLines; {
				activeLineages.add(l);//.getNr());
				String sampleID = tree.getNode(l).getID();
				int sampleState = 0;
				if (states > 1){				
					String[] splits = sampleID.split("_");
					sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
				}
				for (int i = 0; i < states; i++){
					if (i == sampleState){
						linProbs[currPosition] = 1.0;currPosition++;
					}
					else{
						linProbs[currPosition] = 0.0;currPosition++;
					}
				}
			}	
		}
		linProbsLength = newLength;
		// store the node
        if (!saveRamInput.get()) {
        	storeNode(currTreeInterval, currRatesInterval, linProbs, logP, activeLineages, nextTreeEvent, nextRateShift, incomingLines);
        }
    }
    
    private boolean sampleX(int currTreeInterval, int currRatesInterval) {
		int incomingLines = treeIntervals.getLineagesAdded(currTreeInterval);
		linProbsLength = linProbsLength + 1 * states;
		activeLineages.add(incomingLines);
		return treeIntervals.storedLineagesAdded[currTreeInterval] != incomingLines;
		//return treeIntervals.storedLineagesAdded[currRatesInterval] == incomingLines;
    }
    
    private double coalesce(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
    	int coalLines0 = treeIntervals.getLineagesRemoved(currTreeInterval,0);
    	int coalLines1 = treeIntervals.getLineagesRemoved(currTreeInterval,1);
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines0);//.getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines1);//.getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
			System.out.println("daughter lineages at coalescent event not found");
			return Double.NaN;
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
			} else {
				return Double.NEGATIVE_INFINITY;
			}
        }
        
        int lineageToAdd = tree.getNode(coalLines0).getParent().getNr();
        activeLineages.add(lineageToAdd);        
        
        // get the node state probabilities
		DoubleMatrix pVec = new DoubleMatrix();
		pVec.copy(lambda);
		pVec = pVec.div(pVec.sum());
		
		stateProbabilities[tree.getNode(coalLines0).getParent().getNr() - nrSamples] = pVec;
		
		//double[] linProbsNew  = new double[linProbsLength - states];
		
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
		double [] tmp = linProbs;
		linProbs = linProbsNew;	
		linProbsNew = linProbs;
		linProbsLength = linProbsLength - states;
		
		
		//Remove daughter lineages from the line state probs
		if (daughterIndex1>daughterIndex2){
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex1);
			activeLineages.remove(daughterIndex2);			
		} else {
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex2);
			activeLineages.remove(daughterIndex1);			
		}		
     
		if (lambda.min()<0.0){
			System.err.println("Coalescent probability is: " + lambda.min());
			return Double.NEGATIVE_INFINITY;
		}				
		
		// store the node
        if (!saveRamInput.get()) {
        	storeNode(currTreeInterval, currRatesInterval, linProbs, logP + Math.log(lambda.sum()), activeLineages, nextTreeEvent, nextRateShift, lineageToAdd);
        }
		
		if (lambda.sum()==0)
			return Double.NEGATIVE_INFINITY;
		else
			return Math.log(lambda.sum());
    }
     
    private boolean coalesceX(int currTreeInterval, int currRatesInterval) {
    	Integer coalLines0 = treeIntervals.getLineagesRemoved(currTreeInterval,0);
    	Integer coalLines1 = treeIntervals.getLineagesRemoved(currTreeInterval,1);

    	if (!activeLineages.remove(coalLines0) || !activeLineages.remove(coalLines1)) {
			System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
			System.out.println("daughter lineages at coalescent event not found");
			throw new RuntimeException("coalesceX went wrong at 1");    		
    	}
    	activeLineages.add(tree.getNode(coalLines0).getParent().getNr());        

		linProbsLength = linProbsLength - states;
		
		//return treeIntervals.storedLineagesAdded[currRatesInterval] == activeLineages.get(activeLineages.size() - 1);
		return treeIntervals.storedLineagesAdded[currTreeInterval] != activeLineages.get(activeLineages.size() - 1);

//		
//    	final int daughterIndex1 = activeLineages.indexOf(coalLines0);//.getNr());
//		final int daughterIndex2 = activeLineages.indexOf(coalLines1);//.getNr());
//		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
//			System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
//			System.out.println("daughter lineages at coalescent event not found");
//			throw new RuntimeException("coalesceX went wrong at 1");
//		}
//		DoubleMatrix lambda = DoubleMatrix.zeros(states);
//		
//		/*
//		 * Calculate the overall probability for two strains to coalesce 
//		 * independent of the state at which this coalescent event is 
//		 * supposed to happen
//		 */
//        for (int k = 0; k < states; k++) { 
//        	Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];			
//			if (!Double.isNaN(pairCoalRate)){
//				lambda.put(k, pairCoalRate);
//			}else{
//				throw new RuntimeException("coalesceX went wrong at 2");
//			}
//        }
//        
//        activeLineages.add(tree.getNode(coalLines0).getParent().getNr());        
//        
//        // get the node state probabilities
//		DoubleMatrix pVec = new DoubleMatrix();
//		pVec.copy(lambda);
//		pVec = pVec.div(pVec.sum());
//		
//		stateProbabilities[tree.getNode(coalLines0).getParent().getNr() - nrSamples] = pVec;
//		
//		//double[] linProbsNew  = new double[linProbsLength - states];
//		
////		int linCount = 0;		
////		// add all lineages execpt the daughter lineage to the new p array
////		for (int i = 0; i <= nrLineages; i++){
////			if (i != daughterIndex1 && i != daughterIndex2){
////				for (int j = 0; j < states; j++){
////					linProbsNew[linCount*states + j] = linProbs[i*states + j];
////				}
////				linCount++;
////			}
////		}
////		// add the parent lineage
////		for (int j = 0; j < states; j++){
////			linProbsNew[linCount*states + j] = pVec.get(j);
////		}
////		// set p to pnew
////		double [] tmp = linProbs;
////		linProbs = linProbsNew;	
////		linProbsNew = linProbs;
//		linProbsLength = linProbsLength - states;
//		
//		
//		//Remove daughter lineages from the line state probs
//		if (daughterIndex1>daughterIndex2){
//			// remove the daughter lineages from the active lineages
//			activeLineages.remove(daughterIndex1);
//			activeLineages.remove(daughterIndex2);			
//		} else {
//			// remove the daughter lineages from the active lineages
//			activeLineages.remove(daughterIndex2);
//			activeLineages.remove(daughterIndex1);			
//		}
//		
//		//return treeIntervals.storedLineagesAdded[currRatesInterval] == activeLineages.get(activeLineages.size() - 1);
//		return treeIntervals.storedLineagesAdded[currTreeInterval] != activeLineages.get(activeLineages.size() - 1);
	}
    
    public DoubleMatrix getStateProb(int nr){
    	return stateProbabilities[nr - nrSamples];
    }    
    
    public DoubleMatrix getRootState(){
    	return stateProbabilities[stateProbabilities.length-1];
    }
    
    public String getType(){
   		return "state";
    }            
    
    private void storeNode(int storingTreeInterval, int storingRatesInterval, double[] storeLinProbs,
    		double probability, ArrayList<Integer> storeActiveLineages, double nextTreeEvent, double nextRateShift,
    		int addedLineage) {


    	coalRatesInterval[storingTreeInterval] = storingRatesInterval;
    	coalLinProbs[storingTreeInterval] = Arrays.copyOf(storeLinProbs, linProbsLength);
    	coalLogP[storingTreeInterval] = probability;
    	if (storingTreeInterval > 0) {
    		nextTreeEvents[storingTreeInterval - 1] = nextTreeEvent;
    		nextRateShifts[storingTreeInterval - 1] = nextRateShift;
    	}
    	coalActiveLineages[storingTreeInterval].clear();
    	coalActiveLineages[storingTreeInterval].addAll(storeActiveLineages);
    }
        
    private int restoreNode(int restoringInterval){
    	//Log.warning("Restoring " + first + " " + restoringInterval);
    	System.arraycopy(coalLinProbs[restoringInterval], 0, linProbs, 0, coalLinProbs[restoringInterval].length);
    	linProbsLength = coalLinProbs[restoringInterval].length;
    	logP = coalLogP[restoringInterval];    	
    	activeLineages = new ArrayList<>(coalActiveLineages[restoringInterval]);
    	return coalRatesInterval[restoringInterval];

    }
    
    @Override
	public void store(){
    	if (saveRamInput.get()){
    		super.store();
    		return;
    	}
    	// store the intermediate results
    	for (int i = 0; i < coalLinProbs.length; i++) {
    		storeLinProbs[i] = Arrays.copyOf(coalLinProbs[i], coalLinProbs[i].length);
    	}
    	System.arraycopy(coalLogP, 0, storeLogP, 0, coalLogP.length);
    	System.arraycopy(coalRatesInterval, 0, storeRatesInterval, 0, coalRatesInterval.length);
        
    	//System.arraycopy(treeIntervalNrs, 0, storedTreeIntervalNrs, 0, treeIntervalNrs.length);

    	System.arraycopy(nextTreeEvents, 0, storedNextTreeEvents, 0, nextTreeEvents.length);
    	System.arraycopy(nextRateShifts, 0, storedNextRateShifts, 0, nextRateShifts.length);
 
    	//System.arraycopy(lineagesAddded, 0, storedLineagesAddded, 0, lineagesAddded.length);

		for (int i= 0; i < coalActiveLineages.length; i++){
			storeActiveLineages[i].clear();
			storeActiveLineages[i].addAll(coalActiveLineages[i]);
		}
		
//        storedFirst = first;
    	// store the 
    	super.store();
    }
        
    @Override
	public void restore(){
    	if (saveRamInput.get()){
    		super.restore();
    		return;
    	}

    	// store the intermediate results
    	double [] tmp = storeLogP;
    	storeLogP = coalLogP;
    	coalLogP = tmp;
    	
    	double [][] tmp3 = coalLinProbs;
    	coalLinProbs = storeLinProbs;
    	storeLinProbs = tmp3; //Arrays.copyOf(storeLinProbs, storeLinProbs.length);

    	int [] tmp2 = coalRatesInterval;
    	coalRatesInterval = storeRatesInterval;
    	storeRatesInterval = tmp2;
    	
    	//tmp2 = treeIntervalNrs;
    	//treeIntervalNrs = storedTreeIntervalNrs;
    	//storedTreeIntervalNrs = tmp2;
    	
    	//tmp2 = lineagesAddded;
    	//lineagesAddded = storedLineagesAddded;
    	//storedLineagesAddded = lineagesAddded;
    	
    	// coalActiveLineages = new ArrayList<>();
		for (int i= 0; i < storeActiveLineages.length; i++){
			coalActiveLineages[i].clear();
			coalActiveLineages[i].addAll(storeActiveLineages[i]);			
		}    	
    	
		tmp = nextTreeEvents;
		nextTreeEvents = storedNextTreeEvents;
		storedNextTreeEvents = tmp;
		
		tmp = nextRateShifts;
		nextRateShifts = storedNextRateShifts;
		storedNextRateShifts = tmp;

//    	first = storedFirst;    	
        // store the 
    	super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        return ((CalculationNode) dynamics).isDirtyCalculation() || super.requiresRecalculation();
    }

    
}
