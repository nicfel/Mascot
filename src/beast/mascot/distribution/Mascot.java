package beast.mascot.distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math4.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math4.ode.FirstOrderIntegrator;
import org.apache.commons.math4.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math4.ode.nonstiff.HighamHall54Integrator;
import org.jblas.DoubleMatrix;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.IntervalType;
import beast.mascot.dynamics.Dynamics;
import beast.mascot.ode.Euler2ndOrder;
import beast.mascot.ode.MascotODE;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
public class Mascot extends StructuredTreeDistribution {
	
	public Input<Dynamics> dynamicsInput = new Input<>("dynamics", "Input of rates", Input.Validate.REQUIRED);
	public Input<Integer> eulerSteps = new Input<>("steps", "number of steps for euler integration, if not " +
													"specified, dormand prince integration is used");
	public Input<Double> epsilonInput = new Input<>("epsilon", "step size for the RK4 integration",0.001);
	public Input<Double> maxStepInput = new Input<>("maxStep", "step size for the RK4 integration", Double.POSITIVE_INFINITY);
	public Input<Double> stepSizeInput = new Input<>("stepSize", "step size for the RK4 integration");
	public Input<Boolean> saveRamInput = new Input<>("saveRamInput", "doesn't save intermediate steps", false);
	
    
    
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
    private double maxTolerance = 1e-3;            
    private boolean recalculateLogP;    
    
    private boolean useEuler = false;
    
    private double[] times;   
        
    @Override
    public void initAndValidate(){    	
    	treeIntervalsInput.get().calculateIntervals();       
    	stateProbabilities = new DoubleMatrix[treeIntervalsInput.get().getSampleCount()];
        nrSamples = treeIntervalsInput.get().getSampleCount() + 1;    
        states = dynamicsInput.get().getDimension();
        
        times = new double[5];
        
    	int intCount = treeIntervalsInput.get().getIntervalCount();

    	// initialize storing arrays and ArrayLists
    	coalLinProbs = new double[intCount][];
    	coalLogP = new double[intCount];
    	coalRatesInterval = new int[intCount];
    	coalActiveLineages = new ArrayList<>();
    	ArrayList<Integer> emptyList = new ArrayList<>();
    	for (int i = 0; i <= intCount; i++) coalActiveLineages.add(emptyList);
    	
    	if (eulerSteps.get()!=null) useEuler = true;
    	
    	
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
		coalescentRates = dynamicsInput.get().getCoalescentRate(ratesInterval);  
        migrationRates = dynamicsInput.get().getBackwardsMigration(ratesInterval);
        // Time to the next rate shift or event on the tree
        double nextTreeEvent = treeIntervalsInput.get().getInterval(treeInterval);
        double nextRateShift = dynamicsInput.get().getInterval(ratesInterval);
        
        if (!saveRamInput.get()){
	    	if (!first) { 	
	    		// Check if the first interval is already dirty
				if (!dynamicsInput.get().intervalIsDirty(ratesInterval) 
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
			    				coalescentRates = dynamicsInput.get().getCoalescentRate(ratesInterval);  
			    		        migrationRates = dynamicsInput.get().getBackwardsMigration(ratesInterval);
			     	       		nrLineages = activeLineages.size();
			    				break;
			    			}
		            	}
		            	else{
			    			// check if next interval is dirty, if yes stay in that gene interval
			    			if (dynamicsInput.get().intervalIsDirty(ratesInterval+1)){
			    				// use the next species time from the point of the last gene coalescent event
			    				nextRateShift = lastRatesTime;
			    				nextTreeEvent = treeIntervalsInput.get().getInterval(treeInterval);
				    			ratesInterval = lastRatesInterval;
			    				// restore from last (gene) coalescent event
			    				restoreNode(treeInterval);
			    				coalescentRates = dynamicsInput.get().getCoalescentRate(ratesInterval);  
			    		        migrationRates = dynamicsInput.get().getBackwardsMigration(ratesInterval);
			     	       		
			     	       		nrLineages = activeLineages.size();
			    				break;
			    			}
			    			states--;
			    			ratesInterval++; 
			    			nextTreeEvent -= nextRateShift;
			    			nextRateShift = dynamicsInput.get().getInterval(ratesInterval);
		            	}
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
    		    	recalculateLogP = false;
    				System.exit(0);
                	return calculateLogP();
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
	        		Euler2ndOrder euler = new Euler2ndOrder(migrationRates, coalescentRates, nrLineages , coalescentRates.length, epsilonInput.get(), maxStepInput.get());
		        	double[] linProbs_tmp = new double[linProbs.length+1]; 
		        	double[] linProbs_tmpdt = new double[linProbs.length+1]; 
		        	double[] linProbs_tmpddt = new double[linProbs.length+1]; 
		        	double[] linProbs_tmpdddt = new double[linProbs.length+1]; 
		        	
		        	for (int i = 0; i < linProbs.length; i++) linProbs_tmp[i] = linProbs[i];		        	
		        	
		        	linProbs[linProbs.length-1] = 0;
		        	euler.calculateValues(nextEventTime, linProbs_tmp, linProbs_tmpdt, linProbs_tmpddt, linProbs_tmpdddt, times);		        	
	        		
		            for (int i = 0; i < linProbs.length; i++) linProbs[i] = linProbs_tmp[i]; 
		            
		            logP += linProbs_tmp[linProbs_tmp.length-1];
	        	}        	
        	}
       	
        	if (nextTreeEvent <= nextRateShift){
 	        	if (treeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.COALESCENT) {
	        		logP += normalizeLineages();									// normalize all lineages before event		
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		logP += coalesce(treeInterval, ratesInterval);	  				// calculate the likelihood of the coalescent event
	        	}
 	       		
 	       		if (treeIntervalsInput.get().getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	       			
 	       			if (linProbs.length > 0)
 	       				logP += normalizeLineages();								// normalize all lineages before event
	       			nrLineages++;													// sampling event increases the number of lineages by one
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
        		coalescentRates = dynamicsInput.get().getCoalescentRate(ratesInterval);  
                migrationRates = dynamicsInput.get().getBackwardsMigration(ratesInterval);
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = dynamicsInput.get().getInterval(ratesInterval);
        	}
        }while(nextTreeEvent <= Double.POSITIVE_INFINITY);
        first = false;
//        System.exit(0);
        
//        System.out.println(times[1]/(times[0]+times[1]+times[2]+times[3]+times[4]));
//        System.out.println(times[0]/(times[0]+times[1]) + " " + times[2]/(times[2]+times[3]));
        
//        System.out.println((times[0]/(times[0]+times[1]+times[2]+times[3])) + " " + times[1]/(times[0]+times[1]+times[2]+times[3]) +
//        		" " + times[2]/(times[0]+times[1]+times[2]+times[3]) + " " + times[3]/(times[0]+times[1]+times[2]+times[3]));
		return logP;  	
    }   
    
//    private void integrate(double duration){
//    }
//
//    private void ei(double duration, double[] linProbs_for_ode, double[] meanLinProbs){
//    	eulerIntegration(duration, linProbs_for_ode, meanLinProbs);   	
//    }

    
    private void eulerIntegration(double duration, double[] linProbs_for_ode, double[] meanLinProbs){
    	int nrSteps = eulerSteps.get();
    	double integrationWidth = duration/nrSteps;
    	
    	double[] maxRateOut = new double[states];
    	for (int a = 0; a < states; a++)
    		for (int b = 0; b < states; b++)
    			if (a==b)
    				maxRateOut[a] += coalescentRates[a];
    			else
    				maxRateOut[a] += migrationRates[a][b];
    	
    	double max = maxRateOut[0];
    	for (int i = 1; i < states; i++)
    		if (maxRateOut[i]>max)
    			max = maxRateOut[i];
    			
    	double maxStepSize = 1/(1*max);
    	
    	while (integrationWidth > maxStepSize){
    		nrSteps++;
    		integrationWidth = duration/nrSteps;
    	}    
    	
//    	System.out.println(integrationWidth);
    	

    	normalizeLineages();
    	
		for (int intStep = 0; intStep < nrSteps; intStep++){
        	double[] sumStates = new double[states];
        	
        	for (int i = 0; i< nrLineages; i++)
        		for (int j = 0; j<states; j++)
        			sumStates[j] += linProbs[states*i+j];  
        	
        	updateLinProbs(sumStates, linProbs, linProbs_for_ode, integrationWidth);
        	
        	// normalize linProbs_for_ode
        	normalizeLinProbs_for_ode(linProbs_for_ode);
        	getMean(meanLinProbs, linProbs_for_ode);
        	
        	
        	
        	// update logP
        	double[] meanSumStates = new double[states];
        	double[] dTdtStates = new double[states];
        	
        	for (int i = 0; i< nrLineages; i++)
        		for (int j = 0; j<states; j++)
        			meanSumStates[j] += meanLinProbs[states*i+j];  
        	
        	for (int i = 0; i<nrLineages; i++)
        		for (int j = 0; j<states; j++)
        			dTdtStates[j] += meanLinProbs[i*states+j]*(meanSumStates[j] - meanLinProbs[i*states+j]);
        	
    		for (int j = 0; j<states; j++)
    			logP -= integrationWidth*coalescentRates[j]*dTdtStates[j];
    		
            for (int i = 0; i<linProbs.length; i++)
        		linProbs[i] = linProbs_for_ode[i]; 
    	}
    	
    }
    
    private void getMean(double[] meanLinProbs, double[] linProbs_for_ode) {
    	for (int i = 0; i < nrLineages; i++)
    		for (int j = 0; j < states; j++)
    			meanLinProbs[i*states+j] = (linProbs_for_ode[i*states+j] + linProbs[i*states+j])/2;
	}

	private void normalizeLinProbs_for_ode(double[] linProbs_for_ode){
    	for (int i = 0; i < nrLineages; i++){
    		double lineProbs = 0.0;
    		for (int j = 0; j < states; j++)
    			lineProbs += linProbs_for_ode[i*states+j];
    		for (int j = 0; j < states; j++)
    			linProbs_for_ode[i*states+j] = linProbs_for_ode[i*states+j]/lineProbs;
    	}       
    }
    
    private void updateLinProbs(double[] sumStates, double[] linProbs, double[] linProbs_for_ode, double integrationWidth){
//		double[] linProbs_for_ode = new double[linProbs.length];
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<nrLineages; i++){
    		double totCoalRate = getTotCoal(i, sumStates);
    		
   		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			
    			double migrates = getMigration(i, j);
    			
    			// Calculate the Derivate of p:
    			double pDot = migrates + linProbs[states*i+j] *(totCoalRate -2*coalescentRates[j] * (sumStates[j] - linProbs[states*i+j]));
    			double newVal = linProbs[states*i+j] + pDot*integrationWidth;
    			
    			
    			if (newVal < 0.0 || newVal > 1.0){
//    				System.err.println("invalid lin Prob value " + newVal);
    				recalculateLogP = true;
    			}else{
    				linProbs_for_ode[states*i+j] = newVal;
    			}
    		}
    	}   
    }
    
    private double getTotCoal(int i, double[] sumStates){
		double totCoalRate = 0.0;
		for (int j = 0; j<states; j++)
			totCoalRate += 2*coalescentRates[j] * linProbs[states*i+j]* (sumStates[j] - linProbs[states*i+j]);     		
		return totCoalRate;
    }

    private double getMigration(int i, int j){
    	double migrates = 0.0;
		for (int k = 0; k<states; k++){
			if (j != k){				    					
				// the probability of lineage i being in state j is p[i*nr_states +j]
				migrates += linProbs[states*i+k]*migrationRates[k][j] -
						linProbs[states*i+j]*migrationRates[j][k];
			}
		}   	
		return migrates;
    }
    
    private double normalizeLineages(){
    	if (linProbs==null)
    		return 0.0;
    	
    	
    	double interval = 0.0;
    	for (int i = 0; i < nrLineages; i++){
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
		return Math.log(interval/(nrLineages));

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
        if (saveRamInput.get()) return;


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
    	if (saveRamInput.get()){
    		super.store();
    		return;
    	}
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
        // store the 
    	super.restore();
    }

    
}
