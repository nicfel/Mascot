package beast.mascot.distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math4.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math4.ode.FirstOrderIntegrator;
import org.apache.commons.math4.ode.nonstiff.DormandPrince853Integrator;
import org.jblas.DoubleMatrix;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.IntervalType;
import beast.mascot.dynamics.Dynamics;
import beast.mascot.ode.MascotODE;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
public class Mascot extends StructuredTreeDistribution {
	
	public Input<Dynamics> structuredRateIntervals = new Input<>("rates", "Input of rates", Input.Validate.REQUIRED);
//    public Input<RatesAndTraits> ratesAndTraitsInput = new Input<>("ratesAndTraits", "particle model for input");          
	
    
    
	public int samples;
	public int nrSamples;
	public DoubleMatrix[] stateProbabilities;
    
    private int nrLineages;   

    // current rates         
    private double[][] migrationRates;
    private double[] coalescentRates; 	

    
    // Set up for lineage state probabilities
    private ArrayList<Integer> activeLineages;
	private double[] linProbs;
	private int states;
	
    // store the linProbs, multiplicators and logP's at coalescent points in jagged arrays from last time
    private double[][] coalLinProbs;
    private double[] coalLogP;
    private int[] coalRatesInterval;
    private ArrayList<ArrayList<Integer>> coalActiveLineages;
    
    // deep store the things above for MCMC
    private double[][] storeLinProbs;
    private double[] storeLogP;
    private int[] storeRatesInterval;
    private ArrayList<ArrayList<Integer>> storeActiveLineages;
    private boolean storedFirst;

    // check if this is the first calculation
    private boolean first = true;

	
	// maximum integration error tolerance
    private double maxTolerance = 1e-2;            
    private boolean recalculateLogP;    
   
        
    @Override
    public void initAndValidate(){    	
    	treeIntervalsInput.get().calculateIntervals();       
    	stateProbabilities = new DoubleMatrix[treeIntervalsInput.get().getSampleCount()];
        nrSamples = treeIntervalsInput.get().getSampleCount() + 1;    
        states = structuredRateIntervals.get().getDimension();
        
    	int intCount = treeIntervalsInput.get().getIntervalCount();

    	// initialize storing arrays and ArrayLists
    	coalLinProbs = new double[intCount][];
    	coalLogP = new double[intCount];
    	coalRatesInterval = new int[intCount];
    	coalActiveLineages = new ArrayList<>();
    	ArrayList<Integer> emptyList = new ArrayList<>();
    	for (int i = 0; i <= intCount; i++) coalActiveLineages.add(emptyList);
    }
        
    public double calculateLogP() {
    	// newly calculate tree intervals
    	treeIntervalsInput.get().calculateIntervals();
    	// correctly calculate the daughter nodes at coalescent intervals in the case of
    	// bifurcation or in case two nodes are at the same height
    	treeIntervalsInput.get().swap();    	
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages = new ArrayList<Integer>(); 
        logP = 0;
        nrLineages = 0;
        linProbs = new double[0];// initialize the tree and rates interval counter
        int treeInterval = 0, ratesInterval = 0;        
        double nextEventTime = 0.0;
		coalescentRates = structuredRateIntervals.get().getCoalescentRate(ratesInterval);  
        migrationRates = structuredRateIntervals.get().getBackwardsMigration(ratesInterval);
        // Time to the next rate shift or event on the tree
        double nextTreeEvent = treeIntervalsInput.get().getInterval(treeInterval);
        double nextRateShift = structuredRateIntervals.get().getInterval(ratesInterval);
        
    	if (!first) { 	
    		// Check if the first interval is already dirty
			if (!structuredRateIntervals.get().intervalIsDirty(ratesInterval) 
					&& !treeIntervalsInput.get().intervalIsDirty(treeInterval)){ 
				double lastRatesTime = nextTreeEvent;
				int lastRatesInterval = ratesInterval;
				
	    		while(treeInterval < Integer.MAX_VALUE){
	            	if (nextTreeEvent < nextRateShift){
	            		treeInterval++;
		        		// Check if the last interval was reached
		        		if (treeInterval == treeIntervalsInput.get().intervalCount){
		        			logP = coalLogP[coalLogP.length-1];
		        			return logP;
		        		}
		        		ratesInterval -= nextTreeEvent;
		        		lastRatesTime = ratesInterval;
		        		lastRatesInterval = ratesInterval;		        		
		        		nextTreeEvent = treeIntervalsInput.get().getInterval(treeInterval);
		    			// check if this gene interval is dirty
		    			if (treeIntervalsInput.get().intervalIsDirty(treeInterval)){
		    				// if the interval is dirty, restore from last coalescent event
		    				ratesInterval = restoreNode(treeInterval-1);
		    				coalescentRates = structuredRateIntervals.get().getCoalescentRate(ratesInterval);  
		    		        migrationRates = structuredRateIntervals.get().getBackwardsMigration(ratesInterval);
		     	       		nrLineages = activeLineages.size();
		    				break;
		    			}
	            	}
	            	else{
		    			// check if next interval is dirty, if yes stay in that gene interval
		    			if (structuredRateIntervals.get().intervalIsDirty(ratesInterval+1)){
		    				// use the next species time from the point of the last gene coalescent event
		    				nextRateShift = lastRatesTime;
		    				nextTreeEvent = treeIntervalsInput.get().getInterval(treeInterval);
			    			ratesInterval = lastRatesInterval;
		    				// restore from last (gene) coalescent event
		    				restoreNode(treeInterval);
		    				coalescentRates = structuredRateIntervals.get().getCoalescentRate(ratesInterval);  
		    		        migrationRates = structuredRateIntervals.get().getBackwardsMigration(ratesInterval);
		     	       		
		     	       		nrLineages = activeLineages.size();
		    				break;
		    			}
		    			states--;
		    			ratesInterval++; 
		    			nextTreeEvent -= nextRateShift;
		    			nextRateShift = structuredRateIntervals.get().getInterval(ratesInterval);
	            	}
	    		}

			}
		}	        
        // Calculate the likelihood
        do {       
        	nextEventTime = Math.min(nextTreeEvent, nextRateShift); 	 
       	
        	if (nextEventTime > 0) {													// if true, calculate the interval contribution        		
                if(recalculateLogP){
    				System.err.println("ode calculation stuck, reducing tolerance, new tolerance= " + maxTolerance);
    				maxTolerance *=0.9;
    				System.exit(0);
                	return calculateLogP();
                }
	        	double[] linProbs_for_ode = new double[linProbs.length];                
                FirstOrderIntegrator integrator = new DormandPrince853Integrator(1e-32, 1e10, maxTolerance, 1e-100);
                // set the maximal number of evaluations
                integrator.setMaxEvaluations((int) 1e5);
               // set the odes
                FirstOrderDifferentialEquations ode = new MascotODE(migrationRates, coalescentRates, nrLineages , coalescentRates.length);
                
                // integrate	        
                try {
                	integrator.integrate(ode, 0, linProbs, nextEventTime, linProbs_for_ode);
                }catch(Exception e){
                	System.out.println("expection");
                	recalculateLogP = true;    				
                }                
                for (int i = 0; i<linProbs.length; i++)
                		linProbs[i] = linProbs_for_ode[i];  
        	}
        	
        	if (nextTreeEvent <= nextRateShift){
 	        	if (treeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.COALESCENT) {
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		logP += normalizeLineages();									// normalize all lineages before event		
        			logP += coalesce(treeInterval, ratesInterval);	  				// calculate the likelihood of the coalescent event
	        	}
 	       		
 	       		if (treeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	       			
 	       			nrLineages++;													// sampling event increases the number of lineages by one
 	       			if (linProbs.length > 0)
 	       				logP += normalizeLineages();								// normalize all lineages before event
 	       			sample(treeInterval, ratesInterval);							// calculate the likelihood of the sampling event if sampling rate is given
	       		}	
 	       		
 	       		treeInterval++;
        		nextRateShift -= nextTreeEvent;   
        		try{
        			nextTreeEvent = treeIntervalsInput.get().getInterval(treeInterval);
        		}catch(Exception e){
        			break;
        		}
        	}else{
        		ratesInterval++;
        		coalescentRates = structuredRateIntervals.get().getCoalescentRate(ratesInterval);  
                migrationRates = structuredRateIntervals.get().getBackwardsMigration(ratesInterval);
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = structuredRateIntervals.get().getInterval(ratesInterval);
        	}
        }while(nextTreeEvent <= Double.POSITIVE_INFINITY);
        first = false;
        return logP;

  	
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
    
    private void sample(int currTreeInterval, int currRatesInterval) {
		List<Node> incomingLines = treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
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
					linProbsNew[currPosition] = 1.0;currPosition++;
				}
				else{
					linProbsNew[currPosition] = 0.0;currPosition++;
				}
			}
		}	
		linProbs = linProbsNew;
		// store the node
		storeNode(currTreeInterval, currRatesInterval, linProbs, logP, activeLineages);
    }
    
    private double coalesce(int currTreeInterval, int currRatesInterval) {
		List<Node> coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
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
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
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
        	Double pairCoalRate = coalescentRates[k] * 2 * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];			
			if (!Double.isNaN(pairCoalRate)){
				lambda.put(k, pairCoalRate);
			}else{
				return Double.NEGATIVE_INFINITY;
			}
        }
        
        activeLineages.add(coalLines.get(0).getParent().getNr());        
        
        // get the node state probabilities
		DoubleMatrix pVec = new DoubleMatrix();
		pVec.copy(lambda);
		pVec = pVec.div(pVec.sum());
		
		stateProbabilities[coalLines.get(0).getParent().getNr() - nrSamples] = pVec;
		
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
		
		// store the node
		storeNode(currTreeInterval, currRatesInterval, linProbs, logP + Math.log(lambda.sum()), activeLineages);
		
    	return Math.log(lambda.sum());
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
    		double probability, ArrayList<Integer> storeActiveLineages){

    	coalRatesInterval[storingTreeInterval] = storingRatesInterval;
    	coalLinProbs[storingTreeInterval] = Arrays.copyOf(storeLinProbs, storeLinProbs.length);
    	coalLogP[storingTreeInterval] = probability;   	
    	ArrayList<Integer> tmp_activeLins = new ArrayList<>(storeActiveLineages);
    	coalActiveLineages.set(storingTreeInterval, tmp_activeLins);
    }
        
    private int restoreNode(int restoringInterval){   
    	linProbs = Arrays.copyOf(coalLinProbs[restoringInterval], coalLinProbs[restoringInterval].length);
    	logP = coalLogP[restoringInterval];    	
    	activeLineages = new ArrayList<>(coalActiveLineages.get(restoringInterval));
    	return coalRatesInterval[restoringInterval];

    }
    
    @Override
	public void store(){
    	// store the intermediate results
    	storeLinProbs = new double[coalLinProbs.length][];
    	for (int i = 0; i < coalLinProbs.length; i++)
    		storeLinProbs[i] = Arrays.copyOf(coalLinProbs[i], coalLinProbs[i].length);    	
        storeLogP = Arrays.copyOf(coalLogP, coalLogP.length);
        storeRatesInterval = Arrays.copyOf(coalRatesInterval, coalRatesInterval.length);
        
        storeActiveLineages = new ArrayList<>();
		for (int i= 0; i < coalActiveLineages.size(); i++){
			ArrayList<Integer> add = new ArrayList<>();
			for (int j = 0; j < coalActiveLineages.get(i).size(); j++)
				add.add(coalActiveLineages.get(i).get(j));
			storeActiveLineages.add(add);			
		}
		
        storedFirst = first;
//        System.out.println("store........");
    	// store the 
    	super.store();
    }
        
    @Override
	public void restore(){
    	// store the intermediate results
    	coalLinProbs = Arrays.copyOf(storeLinProbs, storeLinProbs.length);
    	coalLogP = Arrays.copyOf(storeLogP, storeLogP.length);   
    	coalRatesInterval = Arrays.copyOf(storeRatesInterval, storeRatesInterval.length);
    	
    	coalActiveLineages = new ArrayList<>();
		for (int i= 0; i < storeActiveLineages.size(); i++){
			ArrayList<Integer> add = new ArrayList<>();
			for (int j = 0; j < storeActiveLineages.get(i).size(); j++)
				add.add(storeActiveLineages.get(i).get(j));
			coalActiveLineages.add(add);			
		}    	
    	
    	first = storedFirst;
    	
//        System.out.println("restore........");
        // store the 
    	super.restore();
    }

    
}