package beast.mascot.distribution;


import java.util.ArrayList;

import org.jblas.DoubleMatrix;

import beast.core.Citation;
import beast.core.Description;
import beast.evolution.tree.TreeInterface;
import beast.mascot.dynamics.Dynamics;
import beast.mascot.ode.*;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
@Citation("Nicola F. Müller, David A. Rasmussen, Tanja Stadler (2017)\n  The Structured Coalescent and its Approximations.\n  Mol Biol Evol 2017 msx186. doi: 10.1093/molbev/msx186")
public class MascotNative extends StructuredTreeDistribution {
	
    
	public int samples;
	public int nrSamples;
	public double[][] stateProbabilities;
    
    private int nrLineages;   

    // current rates         
    private double[] coalescentRates; 	

    
    // Set up for lineage state probabilities
    private ArrayList<Integer> activeLineages;
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
    private int first = 0;

	
	// maximum integration error tolerance
    private double maxTolerance = 1e-3;            
    private boolean recalculateLogP;
	Euler2ndOrderBase euler;
	StructuredTreeIntervals treeIntervals;

	int [] nodeType;

	int [] parents;
	int [] lineagesAdded;
	int [] lineagesRemoved;
	double [] intervals;
	
	int [] storedLineagesAdded;
	int intervalCount;
	
    public MascotNative(Euler2ndOrderBase euler, StructuredTreeIntervals treeIntervals, int [] nodeType, int states) {
    	this.euler = euler;
    	this.treeIntervals = treeIntervals;
    	TreeInterface tree = treeIntervals.treeInput.get();
    	this.nodeType = nodeType;
    	this.states = states;
    	parents = new int[tree.getNodeCount()];

    	treeIntervals.calculateIntervals();       
    	stateProbabilities = new double[treeIntervals.getSampleCount()][states];
        nrSamples = treeIntervals.getSampleCount() + 1;    
                
    	intervalCount = treeIntervals.getIntervalCount();

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
    	//coalActiveLineages = new ArrayList<>();
    	nextTreeEvents = new double[intervalCount];
    	nextRateShifts = new double[intervalCount];
    	storedNextTreeEvents = new double[intervalCount];
    	storedNextRateShifts = new double[intervalCount];
    	//treeIntervalNrs = new int[intCount];
    	//storedTreeIntervalNrs = new int[intCount];
    	//lineagesAddded = new int[intCount];
    	//storedLineagesAddded = new int[intCount];
    	
    	//ArrayList<Integer> emptyList = new ArrayList<>();
    	//for (int i = 0; i <= intCount; i++) coalActiveLineages.add(emptyList);
    	
    	activeLineages = new ArrayList<>();

    	int MAX_SIZE = intervalCount * states;
    	linProbs_for_ode = new double[MAX_SIZE];
    	linProbs_tmp = new double[MAX_SIZE];
    	linProbs = new double[MAX_SIZE];
    	linProbsNew = new double[MAX_SIZE];
    }
    
    double [] linProbs_for_ode;
    double [] linProbs_tmp;
    

    public double calculateLogP(boolean useCache) {
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages.clear();
        logP = 0;
        nrLineages = 0;
        //linProbs = new double[0];// initialize the tree and rates interval counter
        linProbsLength = 0;
        int treeInterval = 0, ratesInterval = 0;        
     	double nextEventTime = 0.0;
        
        // Time to the next rate shift or event on the tree
        double nextTreeEvent = intervals[treeInterval];//treeIntervals.getInterval(treeInterval);
        double nextRateShift = getRateShiftInterval(ratesInterval);
        
        if (first > 0 && useCache) {
        // restore the likelihood to last known good place
    	  int pos0 = -1, pos1 = -1;
    	  do {
        	
        	nextEventTime = nextTreeEvent;

    		// Check if the last interval was reached
    		if (treeInterval == intervalCount){
    			// Log.warning("Restoring to the finish!");
    			logP = coalLogP[coalLogP.length-1];
    			return logP;
    		}
    		boolean isDirty;
    		if (lineagesRemoved[treeInterval*2] >= 0) { // == IntervalType.COALESCENT) {        		
 	           	Integer coalLines0 = lineagesRemoved[treeInterval * 2 + 0];
 	           	Integer coalLines1 = lineagesRemoved[treeInterval * 2 + 1];
 	           	pos0 = activeLineages.indexOf(coalLines0);
 	           	pos1 = activeLineages.indexOf(coalLines1);
 	           	if (pos0 < 0 || pos1 < 0) {
 	           		System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
 	           		System.out.println("daughter lineages at coalescent event not found");
 	           		throw new RuntimeException("coalesceX went wrong at 1");
 	           	}
 	           	if (pos0 > pos1) {
 	           		activeLineages.remove(pos0);
 	           		activeLineages.remove(pos1);
 	           	} else {
 	           		activeLineages.remove(pos1);
 	           		activeLineages.remove(pos0);
 	           	}
 	            int newLineage = parents[coalLines0];
 	            activeLineages.add(newLineage);	 	           
 	            isDirty = storedLineagesAdded[treeInterval] != newLineage;
        	} else { // (treeIntervals.getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	       			
       			int incomingLines = lineagesAdded[treeInterval];
       			activeLineages.add(incomingLines);
       			isDirty = storedLineagesAdded[treeInterval] != incomingLines;
       		}	
       		
       	
    		if (treeIntervals.intervalIsDirty(treeInterval+1) || isDirty) { 
    			if (treeInterval <= 2) {
        			treeInterval = 0;
        			ratesInterval = 0;
        	        activeLineages.clear();
        	        logP = 0;
        	     	nextEventTime = 0.0;
        	        nextTreeEvent = intervals[treeInterval];
        	        nextRateShift = getRateShiftInterval(ratesInterval);        			
        			break; 	    				
    			}
    			
// 	    			Log.warning.print("Restore I ");
	 	        activeLineages.remove(activeLineages.size() - 1);	 	           
        		if (lineagesRemoved[treeInterval*2] > 0) { // == IntervalType.COALESCENT) {        		
	 	           	Integer coalLines0 = lineagesRemoved[treeInterval * 2 + 0];
	 	           	Integer coalLines1 = lineagesRemoved[treeInterval * 2 + 1];
	 	           	if (pos0 > pos1) {
	 	           		activeLineages.add(pos1, coalLines1);
	 	           		activeLineages.add(pos0, coalLines0);
	 	           	} else {
	 	           		activeLineages.add(pos0, coalLines0);
	 	           		activeLineages.add(pos1, coalLines1);
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
    		try{
    			nextTreeEvent = intervals[treeInterval];
    		} catch(Exception e) {
    			break;
    		}
//        	} else {
//        		ratesInterval++;
//        		nextTreeEvent -= nextRateShift;
// 	       		nextRateShift = getRateShiftInterval(ratesInterval);
//        	}
        } while(nextTreeEvent <= Double.POSITIVE_INFINITY);
    }


		coalescentRates = getCoalescentRate(ratesInterval);  
        //migrationRates = dynamics.getBackwardsMigration(ratesInterval);
		//indicators = dynamics.getIndicators(ratesInterval);
		nrLineages = activeLineages.size();
		linProbsLength = nrLineages * states;

        // Calculate the likelihood
        do {       
        	nextEventTime = Math.min(nextTreeEvent, nextRateShift);       	
        	if (nextEventTime > 0) {													// if true, calculate the interval contribution        		
                if (recalculateLogP) {
    				System.err.println("ode calculation stuck, reducing tolerance, new tolerance= " + maxTolerance);
    				maxTolerance *=0.9;
    		    	recalculateLogP = false;
    				System.exit(0);
                	return calculateLogP();
                }
        		logP += doEuler(nextEventTime, ratesInterval);
        	}
       	
        	if (nextTreeEvent <= nextRateShift){
        		if (lineagesRemoved[treeInterval*2] >= 0) { // == IntervalType.COALESCENT) {        		
// 	        	if (treeIntervals.getIntervalType(treeInterval) == IntervalType.COALESCENT) {
// 	        		System.out.print(String.format("%.3f ", nextTreeEvent));
//	        		logP += normalizeLineages(linProbs);									// normalize all lineages before event		
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		logP += coalesce(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);	  				// calculate the likelihood of the coalescent event
	        	} else {
//	       		
// 	       		if (treeIntervals.getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	       			
 	       			//if (linProbsLength > 0)
 	       			//	logP += normalizeLineages(linProbs);								// normalize all lineages before event
	       			nrLineages++;													// sampling event increases the number of lineages by one
 	       			sample(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);							// calculate the likelihood of the sampling event if sampling rate is given
	       		}	
 	       		
 	       		treeInterval++;
        		nextRateShift -= nextTreeEvent;   
        		try {
        			nextTreeEvent = intervals[treeInterval];
        		} catch(Exception e) {
        			break;
        		}
        	} else {
        		ratesInterval++;
        		coalescentRates = getCoalescentRate(ratesInterval);  
                //migrationRates = dynamics.getBackwardsMigration(ratesInterval);
        		//indicators = dynamics.getIndicators(ratesInterval);  
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = getRateShiftInterval(ratesInterval);
        	}
        	if (logP == Double.NEGATIVE_INFINITY) {
        		return logP;
        	}
        } while(nextTreeEvent <= Double.POSITIVE_INFINITY);

        first++;
//        System.err.println(logP);
		return logP;  	
    }

	private double[] getCoalescentRate(int i) {
    	if (i >= coalescentRatesCache.length) {
    		return coalescentRatesCache[coalescentRatesCache.length-1];
    	} else {
			return coalescentRatesCache[i];
    	}
	}

	private double getRateShiftInterval(int i) {
    	if (i >= nextRateShiftCache.length) {
    		return Double.POSITIVE_INFINITY;
    	} else {
			return nextRateShiftCache[i];
    	}
	}

	double [][] coalescentRatesCache;
	double [][] migrationRatesCache;
	int [][] indicatorsCache;
	double [] nextRateShiftCache;
	
	void setUpDynamics(Dynamics dynamics) {
    	int n = dynamics.getEpochCount();
    	if (coalescentRatesCache == null) {
        	coalescentRatesCache = new double[n][];
        	migrationRatesCache = new double[n][];
        	indicatorsCache = new int[n][];
    	}
    	for (int i = 0; i < n; i++) {
    		coalescentRatesCache[i] = dynamics.getCoalescentRate(i);  
            migrationRatesCache[i] = dynamics.getBackwardsMigration(i);
    		indicatorsCache[i] = dynamics.getIndicators(i);
    	}
    	nextRateShiftCache = dynamics.getIntervals();
    	dynamics.setDynamicsKnown();
		euler.setUpDynamics(coalescentRatesCache, migrationRatesCache, indicatorsCache, nextRateShiftCache);
	}

	double [] storedMigrationRates = new double[0];
    double [] storedCoalescentRates = new double[0];
    int storedNrLineages = -1;
    
	private double doEuler(double nextEventTime, int ratesInterval) {
		//for (int i = 0; i < linProbs.length; i++) linProbs_tmp[i] = linProbs[i];
		if (linProbs_tmp.length != linProbsLength + 1) {
			linProbs_tmp= new double[linProbsLength + 1];
		}
		System.arraycopy(linProbs,0,linProbs_tmp,0,linProbsLength);
		linProbs_tmp[linProbsLength] = 0;

		linProbs[linProbsLength-1] = 0;
		

//		if (dynamics.hasIndicators) {
//			euler.initWithIndicators(migrationRates, indicators, coalescentRates, nrLineages);
//			euler.calculateValues(nextEventTime, linProbs_tmp, linProbsLength + 1);
//		} else {
			euler.initAndcalculateValues(ratesInterval, nrLineages, nextEventTime, linProbs_tmp, linProbsLength + 1);
//		}
		
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
    
//    private double normalizeLineages(double [] linProbs){
//    	if (linProbs==null)
//    		return 0.0;
//    	
//    	
//    	double interval = 0.0;
//    	for (int i = 0; i < nrLineages; i++) {
//    		double lineProbs = 0.0;
//    		int u = i * states;
//    		for (int j = 0; j < states; j++) {
//    			if (linProbs[u]>=0.0){
//    				lineProbs += linProbs[u];
//    			} else {
//    				// try recalculation after lowering the tolerance
//    				System.out.println(linProbs[u]);
//    				recalculateLogP = true;
//    				return Math.log(1.0);
//    			}
//    			u++;
//    		}
//			if (lineProbs==0.0) {
//				return Double.NEGATIVE_INFINITY;
//			}
//			u = i * states;
//    		for (int j = 0; j < states; j++) {
//    			linProbs[u] = linProbs[u]/lineProbs;
//    			u++;
//    		}    		
//    		interval += lineProbs;
//    	}	
//		// return mean P_t(T)
//		return Math.log(interval/(nrLineages));
//
//    }
//    
    private void sample(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
		int incomingLines = lineagesAdded[currTreeInterval];
		int newLength = linProbsLength + 1 * states;
		
		int currPosition = linProbsLength;
		
		/*
		 * If there is no trait given as Input, the model will simply assume that
		 * the last value of the taxon name, the last value after a _, is an integer
		 * that gives the type of that taxon
		 */
		activeLineages.add(incomingLines);//.getNr());
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
       	storeNode(currTreeInterval, currRatesInterval, linProbs, logP, activeLineages, nextTreeEvent, nextRateShift, incomingLines);
    }
          
    private double coalesce(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
    	int coalLines0 = lineagesRemoved[currTreeInterval * 2 + 0];
    	int coalLines1 = lineagesRemoved[currTreeInterval * 2 + 1];
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines0);//.getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines1);//.getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
			System.out.println("daughter lineages at coalescent event not found");
			return Double.NaN;
		}

		double [] lambda = stateProbabilities[parents[coalLines0] - nrSamples];
		
		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */
        for (int k = 0; k < states; k++) { 
        	Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];			
			if (!Double.isNaN(pairCoalRate)){
				lambda[k] = pairCoalRate;
			} else {
				return Double.NEGATIVE_INFINITY;
			}
        }
        
        // get the node state probabilities
        double sum = 0;
        for (double d : lambda) {
       	    sum += d;
        }
        for (int i = 0; i < states; i++) {
       	    lambda[i] /= sum;
        }

        int lineageToAdd = parents[coalLines0];
        activeLineages.add(lineageToAdd);        
        
		
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
			linProbsNew[linCount*states + j] = lambda[j]; // pVec.get(j);
		}
		// set p to pnew
		//double [] tmp = linProbs;
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
		
		double min = Double.POSITIVE_INFINITY;
		for (double d : lambda) {
			if (d < min) {min = d;}
		}
		if (min < 0.0){
			System.err.println("Coalescent probability is: " + min);
			return Double.NEGATIVE_INFINITY;
		}				
	
		// store the node
		storeNode(currTreeInterval, currRatesInterval, linProbs, logP + Math.log(sum), activeLineages, nextTreeEvent, nextRateShift, lineageToAdd);
	
		if (sum==0)
			return Double.NEGATIVE_INFINITY;
		else
			return Math.log(sum);
    }
     
  
    public DoubleMatrix getStateProb(int nr){
    	DoubleMatrix p = new DoubleMatrix(states);
    	for (int i = 0; i < states; i++) {
    		p.put(i, stateProbabilities[nr - nrSamples][i]);
    	}
    	return p;
    }    
    
    public DoubleMatrix getRootState(){
    	DoubleMatrix p = new DoubleMatrix(states);
    	for (int i = 0; i < states; i++) {
    		p.put(i, stateProbabilities[stateProbabilities.length-1][i]);
    	}
    	return p;
    }
    
    public String getType(){
   		return "state";
    }            
    
    private void storeNode(int storingTreeInterval, int storingRatesInterval, double[] storeLinProbs,
    		double probability, ArrayList<Integer> storeActiveLineages, double nextTreeEvent, double nextRateShift,
    		int addedLineage) {
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
    	System.arraycopy(coalLogP, 0, storeLogP, 0, coalLogP.length);
    	System.arraycopy(coalRatesInterval, 0, storeRatesInterval, 0, coalRatesInterval.length);
        

    	System.arraycopy(nextTreeEvents, 0, storedNextTreeEvents, 0, nextTreeEvents.length);
    	System.arraycopy(nextRateShifts, 0, storedNextRateShifts, 0, nextRateShifts.length);
    	
    	System.arraycopy(lineagesAdded, 0, storedLineagesAdded, 0, lineagesAdded.length);
    	
    	super.store();
    }
        
    private void storeLinP() {
    	System.arraycopy(coalLinProbsLengths, 0, storedCoalLinProbsLengths, 0, coalLinProbsLengths.length);
    	System.arraycopy(coalLinProbs, 0, storeLinProbs, 0, coalLinProbsLengths[coalLinProbsLengths.length - 1]);
//    	// store intermediate results
//    	for (int i = 0; i < coalLinProbs.length; i++) {
//    		double [] p = coalLinProbs[i];
//    		double [] q = storeLinProbs[i];
//    		if (p.length == q.length) {
//    			System.arraycopy(p, 0, q, 0, p.length);
//    		} else {
//    			q= Arrays.copyOf(p, p.length);
//    		}
//    	}
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

    
}
