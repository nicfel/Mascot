package beast.mascot.distribution;



import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.State;
import beast.evolution.tree.TreeInterface;
import beast.mascot.dynamics.Dynamics;
import beast.mascot.ode.Euler2ndOrderNative;


/**
 * @author Nicola Felix Mueller
 */

@Description("Native Mascot implementation using JNI")
public class MascotNative extends Distribution {
	
    
	private int nrSamples;
	private double[] stateProbabilities;
    
    private int nrLineages;   

    // current rates         
    private double[] coalescentRates; 	

    
    // Set up for lineage state probabilities
    private int [] activeLineages;
    private int activeLineagesCount;
    
	private double[] linProbs;
	private double[] linProbsNew;
	private int linProbsLength;
	private int states;
	
    // store the linProbs, multiplicators and logP's at coalescent points in jagged arrays from last time
    private double[] coalLinProbs;
    private int [] coalLinProbsLengths;
    private double[] coalLogP;
    private int[] coalRatesInterval;
    
    // deep store the things above for MCMC
    private double[] storeLinProbs;
    private int [] storedCoalLinProbsLengths;
    private double[] storeLogP;
    private int[] storeRatesInterval;
    
	private double [] nextTreeEvents;
	private double [] storedNextTreeEvents;
	private double [] nextRateShifts;
	private double [] storedNextRateShifts;

    // check if this is the first calculation
    private boolean first = true;

	
	// maximum integration error tolerance
    private Euler2ndOrderNative euler;

    private int [] nodeType;

    private int [] parents;
    private int [] lineagesAdded;
    private int [] lineagesRemoved;
    private double [] intervals;
	
    private int [] storedLineagesAdded;
    private int intervalCount;
	
    private double [] linProbs_tmp;

	private double [] coalescentRatesCache;
	private double [] migrationRatesCache;
	private int [][] indicatorsCache;
	private double [] nextRateShiftCache;
	private int rateShiftCount;
	
	private double [] currentCoalescentRates;

    public MascotNative(StructuredTreeIntervals treeIntervals, 
    		int [] nodeType, int states, double epsilon, double max_step
    		) {
    	TreeInterface tree = treeIntervals.treeInput.get();
    	treeIntervals.calculateIntervals();  
    	int sampleCount = treeIntervals.getSampleCount();
    	int nodeCount = tree.getNodeCount();
    	intervalCount = treeIntervals.getIntervalCount();
    	setup(nodeType, states, epsilon, max_step, sampleCount, nodeCount, nodeCount);
    }
    
    void setup(int [] nodeType, int states, double epsilon, double max_step, int sampleCount, int nodeCount, int intervalCount) {
    	this.nodeType = nodeType;
    	this.states = states;
    	parents = new int[nodeCount];
    	
    	stateProbabilities = new double[sampleCount * states];
        nrSamples = sampleCount + 1;    
                

    	// initialize storing arrays and ArrayLists
    	storedLineagesAdded = new int[intervalCount];
    	coalLinProbs = new double[intervalCount * intervalCount];
    	storeLinProbs = new double[intervalCount * intervalCount];
    	coalLinProbsLengths = new int[intervalCount];
    	storedCoalLinProbsLengths = new int[intervalCount];
    	coalLogP = new double[intervalCount];
    	storeLogP = new double[intervalCount];
    	coalRatesInterval = new int[intervalCount];
    	storeRatesInterval = new int[intervalCount];
    	nextTreeEvents = new double[intervalCount];
    	nextRateShifts = new double[intervalCount];
    	storedNextTreeEvents = new double[intervalCount];
    	storedNextRateShifts = new double[intervalCount];
    	
    	activeLineages = new int[nodeCount];
    	activeLineagesCount = 0;

    	int MAX_SIZE = intervalCount * states;
    	linProbs_tmp = new double[MAX_SIZE];
    	linProbs = new double[MAX_SIZE];
    	linProbsNew = new double[MAX_SIZE];

    	euler = new Euler2ndOrderNative();
    	euler.setup(MAX_SIZE, states, epsilon, max_step);
    }
    
    

    public double calculateLogP(boolean dynamicsIsDirty, int firstDirtyInterval, int[] lineagesAdded, int[] lineagesRemoved, double[] intervals, int[] parents) {
    	this.lineagesAdded = lineagesAdded;
    	this.lineagesRemoved = lineagesRemoved;
    	this.intervals = intervals;
    	this.parents = parents;
    	
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineagesCount = 0;
        logP = 0;
        nrLineages = 0;
        linProbsLength = 0;
        int treeInterval = 0, ratesInterval = 0;        
     	double nextEventTime = 0.0;
        
        // Time to the next rate shift or event on the tree
        double nextTreeEvent = intervals[treeInterval];//treeIntervals.getInterval(treeInterval);
        double nextRateShift = getRateShiftInterval(ratesInterval);
        
        if (!first && !dynamicsIsDirty && firstDirtyInterval > 2) {
        // restore the likelihood to last known good place
    	  int pos0 = -1, pos1 = -1;
    	  do {
        	
        	nextEventTime = nextTreeEvent;

    		// Check if the last interval was reached
    		if (treeInterval == intervalCount){
    			logP = coalLogP[intervalCount - 1];
    			return logP;
    		}
    		
    		boolean isDirty;
    		if (lineagesRemoved[treeInterval*2] >= 0) { // == IntervalType.COALESCENT) {        		
 	           	int coalLines0 = lineagesRemoved[treeInterval * 2 + 0];
 	           	int coalLines1 = lineagesRemoved[treeInterval * 2 + 1];
 	           	pos0 = indexOf(coalLines0);
 	           	pos1 = indexOf(coalLines1);
 	           	if (pos0 < 0 || pos1 < 0) {
 	           		System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
 	           		System.out.println("daughter lineages at coalescent event not found");
 	           		throw new RuntimeException("coalesceX went wrong at 1");
 	           	}
 	           	if (pos0 > pos1) {
 	           		removeActiveLineageAt(pos0);
 	           		removeActiveLineageAt(pos1);
 	           	} else {
 	           		removeActiveLineageAt(pos1);
 	           		removeActiveLineageAt(pos0);
 	           	}
 	            int newLineage = parents[coalLines0];
 	            addActiveLineages(newLineage);	 	           
 	            isDirty = storedLineagesAdded[treeInterval] != newLineage;
        	} else { // == IntervalType.SAMPLE) { 	       			
       			int incomingLines = lineagesAdded[treeInterval];
       			addActiveLineages(incomingLines);
       			isDirty = storedLineagesAdded[treeInterval] != incomingLines;
       		}	
       		
       	
    		if (isDirty || treeInterval+1 == firstDirtyInterval) { 
    			if (treeInterval <= 2) {
        			treeInterval = 0;
        			ratesInterval = 0;
        	        activeLineagesCount = 0;
        	        logP = 0;
        	     	nextEventTime = 0.0;
        	        nextTreeEvent = intervals[treeInterval];
        	        nextRateShift = getRateShiftInterval(ratesInterval);        			
        			break; 	    				
    			}
    			

	 	        activeLineagesCount--;
        		if (lineagesRemoved[treeInterval*2] > 0) { // == IntervalType.COALESCENT) {        		
	 	           	int coalLines0 = lineagesRemoved[treeInterval * 2 + 0];
	 	           	int coalLines1 = lineagesRemoved[treeInterval * 2 + 1];
	 	           	if (pos0 > pos1) {
	 	           		addActiveLineages(pos1, coalLines1);
	 	           		addActiveLineages(pos0, coalLines0);
	 	           	} else {
	 	           		addActiveLineages(pos0, coalLines0);
	 	           		addActiveLineages(pos1, coalLines1);
	 	           	}
	       		}	
    			
    			treeInterval++;
    			ratesInterval = restoreNode(treeInterval-2);
    			nextTreeEvent = nextTreeEvents[treeInterval-1];
    			nextRateShift = nextRateShifts[treeInterval-1];
    			treeInterval--;
    			break;
    		}
       		treeInterval++;
    		nextRateShift -= nextTreeEvent;
    		if (treeInterval == intervalCount) {
    			break;
    		}
    		nextTreeEvent = intervals[treeInterval];

        } while (nextTreeEvent <= Double.POSITIVE_INFINITY);
    }


		coalescentRates = getCoalescentRate(ratesInterval);  
		nrLineages = activeLineagesCount;
		linProbsLength = nrLineages * states;

        // Calculate the likelihood
        do {       
        	nextEventTime = Math.min(nextTreeEvent, nextRateShift);       	
        	if (nextEventTime > 0) {													// if true, calculate the interval contribution        		
        		logP += doEuler(nextEventTime, ratesInterval);
        	}
       	
        	if (nextTreeEvent <= nextRateShift){
        		if (lineagesRemoved[treeInterval*2] >= 0) { // == IntervalType.COALESCENT) {        		
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		logP += coalesce(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);	  				// calculate the likelihood of the coalescent event
	        	} else {
	       			nrLineages++;													// sampling event increases the number of lineages by one
 	       			sample(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);							// calculate the likelihood of the sampling event if sampling rate is given
	       		}	
 	       		
 	       		treeInterval++;
 	       		if (treeInterval == intervalCount) {
 	       			break;
 	       		}
        		nextRateShift -= nextTreeEvent;   
       			nextTreeEvent = intervals[treeInterval];
        	} else {
        		ratesInterval++;
        		coalescentRates = getCoalescentRate(ratesInterval);  
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = getRateShiftInterval(ratesInterval);
        	}
        	if (logP == Double.NEGATIVE_INFINITY) {
        		return logP;
        	}
        } while (nextTreeEvent <= Double.POSITIVE_INFINITY);

        first = false;
		return logP;  	
    }

	private void addActiveLineages(final int pos1, final int coalLines1) {
		System.arraycopy(activeLineages, pos1, activeLineages, pos1 + 1, activeLineagesCount - pos1);
		activeLineages[pos1] = coalLines1;
		activeLineagesCount++;		
	}

	private void addActiveLineages(int newLineage) {
		activeLineages[activeLineagesCount++] = newLineage;		
	}

	private void removeActiveLineageAt(int pos0) {
		System.arraycopy(activeLineages, pos0 + 1, activeLineages, pos0, activeLineagesCount - pos0);
		activeLineagesCount--;		
	}

	private int indexOf(int value) {
		for (int i = 0; i < activeLineagesCount; i++) {
			if (activeLineages[i] == value) {
				return i;
			}
		}
		return -1;
	}

	private double[] getCoalescentRate(int i) {
    	if (i >= rateShiftCount) {
        	System.arraycopy(coalescentRatesCache, (rateShiftCount-1) * states , currentCoalescentRates, 0, states);
    		return currentCoalescentRates;
    	} else {
        	System.arraycopy(coalescentRatesCache, i * states , currentCoalescentRates, 0, states);
			return currentCoalescentRates;
    	}
	}

	private double getRateShiftInterval(int i) {
    	if (i >= rateShiftCount) {
    		return Double.POSITIVE_INFINITY;
    	} else {
			return nextRateShiftCache[i];
    	}
	}

	
	void setUpDynamics(Dynamics dynamics) {
    	rateShiftCount = dynamics.getEpochCount();
    	if (coalescentRatesCache == null) {
        	coalescentRatesCache = new double[rateShiftCount * states];
        	currentCoalescentRates = new double[states];
        	migrationRatesCache = new double[rateShiftCount * states * states];
        	if (dynamics.hasIndicators) {
        		indicatorsCache = new int[rateShiftCount][];
        	} else {
        		indicatorsCache = new int[1][];
        	}
    	}
    	for (int i = 0; i < rateShiftCount; i++) {
    		System.arraycopy(dynamics.getCoalescentRate(i), 0, coalescentRatesCache, i*states, states);
    		System.arraycopy(dynamics.getBackwardsMigration(i), 0, migrationRatesCache, i*states*states, states*states);
    	}
    	if (dynamics.hasIndicators) {
        	for (int i = 0; i < rateShiftCount; i++) {
        		indicatorsCache = new int[rateShiftCount][];
        	}
    	}
    	nextRateShiftCache = dynamics.getIntervals();
    	dynamics.setDynamicsKnown();
		euler.setUpDynamics(coalescentRatesCache, migrationRatesCache, indicatorsCache, nextRateShiftCache);
	}

    
	private double doEuler(double nextEventTime, int ratesInterval) {
		if (linProbs_tmp.length != linProbsLength + 1) {
			// though less memory efficient, 
			// this reduces jni call time considerably
			linProbs_tmp = new double[linProbsLength + 1];
		}
		System.arraycopy(linProbs,0,linProbs_tmp,0,linProbsLength);
		linProbs_tmp[linProbsLength] = 0;
		linProbs[linProbsLength-1] = 0;
		
		euler.initAndcalculateValues(ratesInterval, nrLineages, nextEventTime, linProbs_tmp, linProbsLength + 1);

		System.arraycopy(linProbs_tmp,0,linProbs,0,linProbsLength);

		return linProbs_tmp[linProbsLength];
	}

	

    private void sample(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
		int incomingLines = lineagesAdded[currTreeInterval];
		int newLength = linProbsLength + 1 * states;
		
		int currPosition = linProbsLength;
		
		/*
		 * If there is no trait given as Input, the model will simply assume that
		 * the last value of the taxon name, the last value after a _, is an integer
		 * that gives the type of that taxon
		 */
		addActiveLineages(incomingLines);//.getNr());
		int sampleState = nodeType[incomingLines];//dynamics.getValue(tree.getNode(l).getID());
						
		for (int i = 0; i < states; i++){
			if (i == sampleState){
				linProbs[currPosition] = 1.0;currPosition++;
			}
			else{
				linProbs[currPosition] = 0.0;currPosition++;
			}
		}

		linProbsLength = newLength;
		// store the node
       	storeNode(currTreeInterval, currRatesInterval, linProbs, logP, nextTreeEvent, nextRateShift);
    }
          
    private double coalesce(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
    	int coalLines0 = lineagesRemoved[currTreeInterval * 2 + 0];
    	int coalLines1 = lineagesRemoved[currTreeInterval * 2 + 1];
		
    	final int daughterIndex1 = indexOf(coalLines0);//.getNr());
		final int daughterIndex2 = indexOf(coalLines1);//.getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
			System.out.println("daughter lineages at coalescent event not found");
			return Double.NaN;
		}

		int offset = (parents[coalLines0] - nrSamples) * states;
		//double [] lambda = stateProbabilities[parents[coalLines0] - nrSamples];
		
		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */
        for (int k = 0; k < states; k++) { 
        	Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];			
			if (!Double.isNaN(pairCoalRate)){
				stateProbabilities[offset + k] = pairCoalRate;
			} else {
				return Double.NEGATIVE_INFINITY;
			}
        }
        
        // get the node state probabilities
        double sum = 0;
        for (int k = 0; k < states; k++) {
       	    sum += stateProbabilities[offset + k];
        }
        for (int i = 0; i < states; i++) {
        	stateProbabilities[offset + i] /= sum;
        }

        int lineageToAdd = parents[coalLines0];
        addActiveLineages(lineageToAdd);        
        
		
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
			linProbsNew[linCount*states + j] = stateProbabilities[offset + j]; // pVec.get(j);
		}
		// set p to pnew
		linProbs = linProbsNew;	
		linProbsNew = linProbs;
		linProbsLength = linProbsLength - states;
		
		
		//Remove daughter lineages from the line state probs
		if (daughterIndex1>daughterIndex2){
			// remove the daughter lineages from the active lineages
			removeActiveLineageAt(daughterIndex1);
			removeActiveLineageAt(daughterIndex2);			
		} else {
			// remove the daughter lineages from the active lineages
			removeActiveLineageAt(daughterIndex2);
			removeActiveLineageAt(daughterIndex1);			
		}
		
		double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < states; i++) {
			if (stateProbabilities[offset + i] < min) {min = stateProbabilities[offset + i];}
		}
		if (min < 0.0){
			System.err.println("Coalescent probability is: " + min);
			return Double.NEGATIVE_INFINITY;
		}				
	
		// store the node
		sum = Math.log(sum);
		storeNode(currTreeInterval, currRatesInterval, linProbs, logP + sum, nextTreeEvent, nextRateShift);
	
		if (sum==0)
			return Double.NEGATIVE_INFINITY;
		else
			return sum;
    }
     
  
    public DoubleMatrix getStateProb(int nr){
    	DoubleMatrix p = new DoubleMatrix(states);
    	for (int i = 0; i < states; i++) {
    		p.put(i, stateProbabilities[(nr - nrSamples) * states + i]);
    	}
    	return p;
    }    
    
    public DoubleMatrix getRootState(){
    	DoubleMatrix p = new DoubleMatrix(states);
    	for (int i = 0; i < states; i++) {
    		p.put(i, stateProbabilities[(nrSamples-2) * states + i]);
    	}
    	return p;
    }
    
    public String getType(){
   		return "state";
    }            
    
    private void storeNode(int storingTreeInterval, int storingRatesInterval, double[] storeLinProbs,
    		double probability, double nextTreeEvent, double nextRateShift) {
    	coalRatesInterval[storingTreeInterval] = storingRatesInterval;
    	int offset = 0;
    	if (storingTreeInterval > 0) {
    		offset = coalLinProbsLengths[storingTreeInterval-1];
    	}
    	System.arraycopy(storeLinProbs, 0, coalLinProbs, offset, linProbsLength);
    	coalLinProbsLengths[storingTreeInterval] = offset + linProbsLength;
    	coalLogP[storingTreeInterval] = probability;
    	nextTreeEvents[storingTreeInterval] = nextTreeEvent;
    	nextRateShifts[storingTreeInterval] = nextRateShift;
    }
        
    private int restoreNode(int restoringInterval){
    	//Log.warning("Restoring " + first + " " + restoringInterval);
    	int offset = 0;
    	if (restoringInterval > 0) {
    		offset = coalLinProbsLengths[restoringInterval-1];
    	}
    	linProbsLength = coalLinProbsLengths[restoringInterval] - offset;
    	System.arraycopy(coalLinProbs, offset, linProbs, 0, linProbsLength);

    	logP = coalLogP[restoringInterval];    	
    	return coalRatesInterval[restoringInterval];

    }
    
    @Override
	public void store() {
    	storeLinP();
    	System.arraycopy(coalLogP, 0, storeLogP, 0, intervalCount);
    	System.arraycopy(coalRatesInterval, 0, storeRatesInterval, 0, intervalCount);
        

    	System.arraycopy(nextTreeEvents, 0, storedNextTreeEvents, 0, intervalCount);
    	System.arraycopy(nextRateShifts, 0, storedNextRateShifts, 0, intervalCount);
    	
    	System.arraycopy(lineagesAdded, 0, storedLineagesAdded, 0, intervalCount);
    	
    	super.store();
    }
        
    private void storeLinP() {
    	System.arraycopy(coalLinProbsLengths, 0, storedCoalLinProbsLengths, 0, intervalCount);
    	System.arraycopy(coalLinProbs, 0, storeLinProbs, 0, coalLinProbsLengths[intervalCount - 1]);
	}


	@Override
	public void restore(){    	
    	// restore intermediate results
    	double [] tmp = storeLogP;
    	storeLogP = coalLogP;
    	coalLogP = tmp;
    	
    	tmp = coalLinProbs;
    	coalLinProbs = storeLinProbs;
    	storeLinProbs = tmp;
    	
    	int [] tmp2 = coalLinProbsLengths;
    	coalLinProbsLengths = storedCoalLinProbsLengths;
    	storedCoalLinProbsLengths = tmp2;

    	tmp2 = coalRatesInterval;
    	coalRatesInterval = storeRatesInterval;
    	storeRatesInterval = tmp2;
    	
		tmp = nextTreeEvents;
		nextTreeEvents = storedNextTreeEvents;
		storedNextTreeEvents = tmp;
		
		tmp = nextRateShifts;
		nextRateShifts = storedNextRateShifts;
		storedNextRateShifts = tmp;

//      not here: lineagesAdded points to TreeIntervals 		
//		tmp2 = lineagesAdded;
//		lineagesAdded = storedLineagesAdded;
//		storedLineagesAdded = tmp2;
		
    	super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}

    
}
