package mascot.distribution;

import beast.base.core.*;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.util.Randomizer;
import cern.colt.Arrays;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
@Citation("Nicola F. Müller, David A. Rasmussen, Tanja Stadler (2017)\n  The Structured Coalescent and its Approximations.\n  Mol Biol Evol 2017 msx186. doi: 10.1093/molbev/msx186")
public class MappedMascot extends Mascot implements Loggable {

	public Input<Double> maxIntegrationStepMappingInput = new Input<>("maxIntegrationStepMapping",
			"maxIntegrationStepMappingas percentage of tree height", 0.001);
	public Input<Double> maxDiffInput = new Input<>("maxDiff",
			"maximial difference to be considered the same intermediate result", 1e-14);


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
	public Input<Boolean> internalNodesOnlyInput = new Input<Boolean>("internalNodesOnly",
			"If true, migration events are not logged default false",
			false);


	Map<Integer, List<double[]>> intermediateStateProbs;
	Map<Integer, List<Double>> intermediateTimes;

	protected DecimalFormat df;
	protected boolean someMetaDataNeedsLogging;
	protected boolean substitutions = false;
	
	private enum EventType {
		SAMPLE, COALESCENT, RATESHIFT
	}

	public Tree mappedTree;

	List<Integer> activeStates;
	List<List<Double>> usedTimes; // keeps track of the timings used for the intermediate results
	List<EventType> usedTypes; // keeps track of the timings used for the intermediate results	
	double[] migrationRates;
	
	long lastLog=-1;

	@Override
	public void initAndValidate() {
    	super.initAndValidate();
    	

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

	public void calcForLogging(long sample) {
		if (lastLog!=sample) {
			try {
				calculateLogP();
			}catch (Exception e) {
				System.out.println(e);
				System.out.println(calculateLogP());
			}
			lastLog=sample;
		}
	}

	@Override
	public double calculateLogP() {
		
//		System.out.println(treeIntervals.treeInput.get());
		// newly calculate tree intervals (already done by swap() below)
		treeIntervals.calculateIntervals();
		// correctly calculate the daughter nodes at coalescent intervals in the case of
		// bifurcation or in case two nodes are at the same height
		treeIntervals.swap();
		
		// initialize a new Mapped tree that inclused single child nodes for migration events		
		mappedTree = new Tree(tree.getRoot().copy());
		mappedTree.getRoot().sort();
		
		// Maps that keep track of the state probabilities and times that fall onto an edge
		intermediateStateProbs = new HashMap<>();
		intermediateTimes = new HashMap<>();
		// reset the usedTimes
		usedTimes = new ArrayList<>();
		usedTypes = new ArrayList<>();

		// The maximum step size for storing intermediate results
		double maxStepSize = treeIntervals.treeInput.get().getRoot().getHeight() * maxIntegrationStepMappingInput.get();
		
		// Set up ArrayLists for the indices of active lineages and the lineage state
		// probabilities
		activeLineages.clear();
		logP = 0;
		nrLineages = 0;
		// linProbs = new double[0];// initialize the tree and rates interval counter
		linProbsLength = 0;
		int treeInterval = 0, ratesInterval = 0, anyInterval = 0;
		double nextEventTime = 0.0;

		// Time to the next rate shift or event on the tree
		double nextTreeEvent = treeIntervals.getInterval(treeInterval);
		double nextRateShift = dynamics.getInterval(ratesInterval);

		if (debug) {
			Log.info.println("##" + treeIntervals.firstDirtyInterval);
			Log.info.println("##" + Arrays.toString(treeIntervals.lineagesAdded));
			Log.info.println("##" + Arrays.toString(treeIntervals.lineagesRemoved));
			Log.info.println("##" + Arrays.toString(treeIntervals.intervals));
		}
		
//		if (first == 0 || !dynamics.areDynamicsKnown()) {
			setUpDynamics();
//		}

		coalescentRates = dynamics.getCoalescentRate(ratesInterval);
		migrationRates = dynamics.getBackwardsMigration(ratesInterval);
		// indicators = dynamics.getIndicators(ratesInterval);
		nrLineages = activeLineages.size();
		linProbsLength = nrLineages * states;

		double currTime = 0;

		double lastRateShift = currTime;
		
		usedTimes.add(new ArrayList<>());
		
		double maxFloatError = treeIntervals.treeInput.get().getRoot().getHeight()*maxDiffInput.get();


		// Calculate the likelihood
		do {
			nextEventTime = Math.min(nextTreeEvent, nextRateShift);
			if (nextEventTime > maxFloatError) { // if true, calculate the interval contribution
				if (recalculateLogP) {
					System.err.println("ode calculation stuck, reducing tolerance, new tolerance= " + maxTolerance);
					maxTolerance *= 0.9;
					recalculateLogP = false;
					System.exit(0);
					return calculateLogP();
				}
				
				// calculate the intermediate state probabilities until the next event, either a coal, sample or rate shift
				if (nextEventTime < maxStepSize) {
					logP += doEuler(nextEventTime, ratesInterval);
					currTime += nextEventTime;
					// add the current time to the used times
					usedTimes.get(anyInterval).add(currTime);
					// store the intermediate results
					storeIntermediateResults(currTime);
				} else {
					int nrIntermediates = (int) (nextEventTime / maxStepSize);
					double oldCurrTime = currTime;
					for (int i = 0; i < (nrIntermediates + 1); i++) {
						logP += doEuler(nextEventTime / (nrIntermediates + 1), ratesInterval);
						currTime += nextEventTime / (nrIntermediates + 1);
						if (i == nrIntermediates)
							currTime = oldCurrTime + nextEventTime;

						// add the current time to the used times
						usedTimes.get(anyInterval).add(currTime);
						// store the intermediate results
						storeIntermediateResults(currTime);
					}
				}
			}

			usedTimes.add(new ArrayList<>());
			anyInterval++;

			if (nextTreeEvent <= nextRateShift) {
				if (treeIntervals.getIntervalType(treeInterval) == IntervalType.COALESCENT) {
					nrLineages--; // coalescent event reduces the number of lineages by one
					logP += coalesce(treeInterval, ratesInterval, nextTreeEvent, nextRateShift, currTime); // calculate
					usedTimes.get(anyInterval).add(currTime);
					usedTypes.add(EventType.COALESCENT);
																											// the
					// likelihood of the
					// coalescent event
				}

				if (treeIntervals.getIntervalType(treeInterval) == IntervalType.SAMPLE) {
					// if (linProbsLength > 0)
					// logP += normalizeLineages(linProbs); // normalize all lineages before event
					nrLineages++; // sampling event increases the number of lineages by one
					sample(treeInterval, ratesInterval, nextTreeEvent, nextRateShift, currTime); // calculate the
					usedTimes.get(anyInterval).add(currTime);
					usedTypes.add(EventType.SAMPLE);

					// the sampling event if
					// sampling rate is given
				}

				treeInterval++;
				nextRateShift -= nextTreeEvent;
				try {
					nextTreeEvent = treeIntervals.getInterval(treeInterval);
				} catch (Exception e) {
					break;
				}
			} else {
				ratesInterval++;
				coalescentRates = dynamics.getCoalescentRate(ratesInterval);
				migrationRates = dynamics.getBackwardsMigration(ratesInterval);
				// indicators = dynamics.getIndicators(ratesInterval);
				nextTreeEvent -= nextRateShift;
				nextRateShift = dynamics.getInterval(ratesInterval);
				lastRateShift = currTime;
				usedTimes.get(anyInterval).add(currTime);
				usedTypes.add(EventType.RATESHIFT);
			}
			if (logP == Double.NEGATIVE_INFINITY) {
				return logP;
			}
			if (debug) {
				Log.info(treeInterval + " " + ratesInterval + " " + logP);
			}
		} while (nextTreeEvent <= Double.POSITIVE_INFINITY);

		first++;
		
		resample(treeInterval, ratesInterval, lastRateShift);
		
		return logP;
	}
	
	@Override
	protected void setUpDynamics() {
    	int n = dynamics.getEpochCount();
    	double [][] coalescentRates = new double[n][];
    	double [][] migrationRates = new double[n][];
    	int [][] indicators = new int[n][];
    	double [] nextRateShift = dynamics.getIntervals();
    	for (int i = 0; i < n; i++) {
    		coalescentRates[i] = dynamics.getCoalescentRate(i);  
            migrationRates[i] = dynamics.getBackwardsMigration(i);
    		indicators[i] = dynamics.getIndicators(i);
    	}
//    	dynamics.setDynamicsKnown();
		euler.setUpDynamics(coalescentRates, migrationRates, indicators, nextRateShift);
	}

	protected double coalesce(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift,
			double currTime) {
		double logP = super.coalesce(currTreeInterval, currRatesInterval, nextTreeEvent, nextRateShift);
		int coalLines0 = treeIntervals.getLineagesRemoved(currTreeInterval, 0);
		int lineageToAdd = tree.getNode(coalLines0).getParent().getNr();
		addNewLineage(lineageToAdd, currTime);
		return logP;
	}

	protected void sample(int treeInterval, int ratesInterval, double nextTreeEvent, double nextRateShift,
			double currTime) {
		super.sample(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);
		int sampLin = treeIntervals.getLineagesAdded(treeInterval);
		addNewLineage(sampLin, currTime);
	}

	private void addNewLineage(int nr, double time) {
		intermediateStateProbs.put(nr, new ArrayList<>());
		intermediateTimes.put(nr, new ArrayList<>());

		final int daughterIndex1 = activeLineages.indexOf(nr);// .getNr());

		double[] probs = new double[states];

		for (int i = 0; i < states; i++)
			probs[i] = linProbs[daughterIndex1 * states + i];
		
		intermediateStateProbs.get(nr).add(probs);
		intermediateTimes.get(nr).add(time);
	}

	
	private void resample(int treeInterval, int ratesInterval, double lastRateShift) {
		treeInterval--;
		// start by resampling the root
		int rootNr = treeIntervals.getLineagesAdded(treeInterval);
		// the number of lineages at the root is always 1
		nrLineages = 1;
		// the active lineages are the ones that are present at the root
		activeStates = new ArrayList<>();
		// sample rootState from the state probabilities at the root
		int rootState = Randomizer.randomChoicePDF(intermediateStateProbs.get(rootNr).get(0));
		// add the rootState to the active states
		activeStates.add(rootState);
		// keep track of the active lineages
		activeLineages.clear();
		activeLineages.add(rootNr);
        linProbsLength = 0;
        // index to keep track of the used times, start with the last one

		mappedTree.getRoot().setMetaData("location", rootState);
		coalesceDown(treeInterval);

		int anyInterval = usedTimes.size() - 2; // start with the series of times saved in the last interval before the coalescent event
		double currTime = usedTimes.get(anyInterval+1).get(0); // get the time of the root as the time of the coalescent interval of the last usedTimes
		double maxFloatError = currTime*maxDiffInput.get();
		
		// get the time until the next tree event and the next rate shift
		double nextTreeEvent = treeIntervals.getInterval(treeInterval);		
		double nextRateShift = currTime - lastRateShift;
				
		double nextEventTime;
		// Calculate the likelihood
		do {
			nextEventTime = Math.min(nextTreeEvent, nextRateShift);
			if (nextEventTime > maxFloatError) { // if true, calculate the interval contribution
				// sample the migration events that occurred between now and the next event time (sample, coal, or rate shift)			
				sampleMigrationEvents(usedTimes.get(anyInterval));
				currTime -= nextEventTime;
			}
			anyInterval--;
			if (usedTypes.get(anyInterval) != EventType.RATESHIFT) {
				if (usedTypes.get(anyInterval) == EventType.COALESCENT) {
					nrLineages++; // coalescent event reduces the number of lineages by one
					coalesceDown(treeInterval - 1); // calculate the
				}

				if (usedTypes.get(anyInterval) == EventType.SAMPLE) {
					// if (linProbsLength > 0)
					// logP += normalizeLineages(linProbs); // normalize all lineages before event
					nrLineages--; // sampling event increases the number of lineages by one
					sampleDown(treeInterval - 1); // calculate the likelihood of
				}

				treeInterval--;
				nextRateShift -= nextTreeEvent;
				try {
					nextTreeEvent = treeIntervals.getInterval(treeInterval);
				} catch (Exception e) {
					break;
				}
			} else {
				nextTreeEvent -= nextRateShift;
				ratesInterval--;
				if (ratesInterval == 0) {
					nextRateShift = Double.POSITIVE_INFINITY;
				}else {
					nextRateShift = dynamics.getInterval(ratesInterval);
					coalescentRates = dynamics.getCoalescentRate(ratesInterval - 1);
					migrationRates = dynamics.getBackwardsMigration(ratesInterval - 1);
				}
				// indicators = dynamics.getIndicators(ratesInterval);
			}
		} while (treeInterval > 0);
	}
	

	private void sampleMigrationEvents(List<Double> timesInInterval) {
		for (int i = 0; i < activeLineages.size(); i++) {
			sampleMigrationEventsLineage(activeLineages.get(i), i, timesInInterval);
		}
	}

	private void sampleMigrationEventsLineage(Integer nodeNr, int index, List<Double> timesInInterval) {
		double K = -Math.log(Randomizer.nextDouble());
		double I = 0.0;
		double startTime = timesInInterval.get(timesInInterval.size() - 1);
		double endTime = timesInInterval.get(0);
		double currentTime = startTime;

		int currTimeInterval = intermediateTimes.get(nodeNr).indexOf(startTime);
		if (currTimeInterval == -1) {
			throw new IllegalArgumentException("timing not found");
		}

		double[] prob_start = intermediateStateProbs.get(nodeNr).get(currTimeInterval);

		while (currentTime > (endTime+maxDiffInput.get())) {
			
			double[] prob_end = intermediateStateProbs.get(nodeNr).get(currTimeInterval - 1);

			int currState = activeStates.get(index);
			
			double[] integral_state = new double[states];
			
			double dt = currentTime - intermediateTimes.get(nodeNr).get(currTimeInterval - 1);
						
			double sumInt = 0;
			for (int i = 0; i < states; i++) {
				if (i != currState) {
					double rates_start = migrationRates[i * states + currState] * 
							prob_start[i] / prob_start[currState];
					double rates_end = migrationRates[i * states + currState] * 
							prob_end[i] / prob_end[currState];
					
					
					if (prob_start[currState] <= 0 || prob_end[currState] <= 0) {
						integral_state[i] = Double.POSITIVE_INFINITY;
					}else if (prob_end[i] <= 0 || prob_start[i] <= 0){
						integral_state[i] = 0;
					}else{
						integral_state[i] = 0.5 * (rates_end + rates_start) * dt;
					}

					sumInt += integral_state[i];
				}
			}

			if ((I + sumInt) > K) {
				// approximate the height of the new event
				double intermediatePoint = Math.max(0.01, (K - I) / sumInt);
				
				currentTime = currentTime
						- dt * (intermediatePoint);
				
				// update stateprobs for this point
				for (int i = 0; i < states; i++)
					prob_start[i] = (intermediatePoint) * prob_start[i] + (1-intermediatePoint) * prob_end[i];
				
				int newState = -1;
				boolean hasInf = false;
				for (int i = 0; i < integral_state.length; i ++) {
					if (integral_state[i] == Double.POSITIVE_INFINITY && prob_end[i] > 0) {
						newState = i;
						hasInf = true;
					}
				}
				
				// sample migration event
				if (!hasInf)
					newState = Randomizer.randomChoicePDF(integral_state);
				
				if (newState == currState)
					System.exit(0);
				
			
				// ad migration event
				Node n = mappedTree.getNode(nodeNr);								
				Node p = n.getParent();
				
				Node migNode = new Node();
				migNode.setMetaData("location", currState);
				migNode.setHeight(currentTime);
				
								
				migNode.setParent(p);
				migNode.addChild(n);


				p.removeChild(n);
				p.addChild(migNode);
				

				n.setParent(migNode);
				
				activeStates.set(index, newState);

				// reset "timers"
				K = -Math.log(Randomizer.nextDouble());
				I = 0.0;
			} else {
				I += sumInt;
				currTimeInterval--;
				if (currTimeInterval == 0)
					break;

				prob_start = intermediateStateProbs.get(nodeNr).get(currTimeInterval);
				currentTime = intermediateTimes.get(nodeNr).get(currTimeInterval);
			}
		}
	}

	private void coalesceDown(int currTreeInterval) {
		int coalLines0 = treeIntervals.getLineagesRemoved(currTreeInterval, 0);
		int coalLines1 = treeIntervals.getLineagesRemoved(currTreeInterval, 1);

		int lineageToAdd = tree.getNode(coalLines0).getParent().getNr();
		
		

		int currState = activeStates.get(activeLineages.indexOf(lineageToAdd));
		mappedTree.getNode(lineageToAdd).setMetaData("location", currState);
		activeStates.remove(activeLineages.indexOf(lineageToAdd));
		activeLineages.remove(activeLineages.indexOf(lineageToAdd));


		activeLineages.add(coalLines0);
		activeLineages.add(coalLines1);
		activeStates.add(currState);
		activeStates.add(currState);

	}
	

	private void sampleDown(int currTreeInterval) {
		
		int incomingLines = treeIntervals.getLineagesAdded(currTreeInterval);
		mappedTree.getNode(incomingLines).setMetaData("location", activeStates.get(activeLineages.indexOf(incomingLines)));
		mappedTree.getNode(incomingLines).getMetaData("location");

		activeStates.remove(activeLineages.indexOf(incomingLines));
		activeLineages.remove(activeLineages.indexOf(incomingLines));
	}

	private void storeIntermediateResults(double time) {
		for (Integer nodeNr : activeLineages) {
			final int daughterIndex1 = activeLineages.indexOf(nodeNr);// .getNr());
			double[] probs = new double[states];
			for (int i = 0; i < states; i++)
				probs[i] = linProbs[daughterIndex1 * states + i];

			intermediateStateProbs.get(nodeNr).add(probs);
			intermediateTimes.get(nodeNr).add(time);
		}
	}

	@Override
	public void store() {
//		super.store();
	}


	@Override
	public void restore() {
//		super.restore();
	}

	@Override
	public void init(PrintStream out) {
		structuredTreeIntervalsInput.get().treeInput.get().init(out);
	}

	public void log(int nSample, PrintStream out) {
		log((long) nSample, out);
	}

	@Override
	public void log(long nSample, PrintStream out) {
		calcForLogging(nSample);

		List<Function> metadata = parameterInput.get();
		for (int i = 0; i < metadata.size(); i++) {
			if (metadata.get(i) instanceof StateNode) {
				metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
			}
		}
		
		BranchRateModel.Base branchRateModel = clockModelInput.get();
		// write out the log tree with meta data
		out.print("tree STATE_" + nSample + " = ");
		root = mappedTree.getRoot();
		if (internalNodesOnlyInput.get())
			out.print(toNewick(root, metadata, branchRateModel, root.getHeight()));
		else
			out.print(toNewick(root, metadata, branchRateModel));
		out.print(";");
	}

	Node root;

	/**
	 * Appends a double to the given StringBuffer, formatting it using the private
	 * DecimalFormat instance, if the input 'dp' has been given a non-negative
	 * integer, otherwise just uses default formatting.
	 * 
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

	public String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
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
//			System.out.println(node.getMetaData("location"));;
//			System.out.println(node.getHeight());
//			System.out.println(lalala);
			
			buf.append("[&");
			try {
				buf.append(dynamics.typeTraitInput.getName() + "=" + dynamics.getStringStateValue((int) node.getMetaData("location")));
			}catch (Exception e) {
				System.out.println(node.getMetaData("location"));
				System.out.println(e);
			}
			
		
			if (branchRateModel != null) {
				buf.append(",rate=");
				appendDouble(buf, branchRateModel.getRateForBranch(node));
			}
			buf.append(']');


		} else {
			String sampleID = node.getID();
			String[] splits = sampleID.split("_");
			int sampleState;

			if (dynamics.typeTraitInput.get() != null) {
				sampleState = dynamics.getValue(node.getID());
			}

			else {
				sampleState = Integer.parseInt(splits[splits.length - 1]); // samples states (or priors) should
																			// eventually be specified in the XML
			}
			
			if ( node.getMetaData("location")!=null) {
				buf.append("[&");
				buf.append(dynamics.typeTraitInput.getName() + "="
						+ dynamics.getStringStateValue((int) node.getMetaData("location")));
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
	
	public String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel, double lastHeight) {
		StringBuffer buf = new StringBuffer();
		if (node.getLeft() != null) {
			if(node.getChildCount()==2) {
				buf.append("(");
				buf.append(toNewick(node.getLeft(), metadataList, branchRateModel, node.getHeight()));
				if (node.getRight() != null) {
					buf.append(',');
					buf.append(toNewick(node.getRight(), metadataList, branchRateModel, node.getHeight()));
				}
				buf.append(")");
			}else {
				buf.append(toNewick(node.getLeft(), metadataList, branchRateModel, lastHeight));
			}
		} else {
			buf.append(node.getNr() + 1);
		}
		
		if(node.getChildCount()!=1) {

			if (!node.isLeaf()) {
				buf.append("[&");
				try {
					buf.append(dynamics.typeTraitInput.getName() + "=" + dynamics.getStringStateValue((int) node.getMetaData("location")));
				}catch (Exception e){
					System.out.println(e);
					throw new IllegalArgumentException("node not mapped");
					
				}
				
			
				if (branchRateModel != null) {
					buf.append(",rate=");
					appendDouble(buf, branchRateModel.getRateForBranch(node));
				}
				buf.append(']');
	
	
			} else {
				String sampleID = node.getID();
				String[] splits = sampleID.split("_");
				int sampleState;
	
				if (dynamics.typeTraitInput.get() != null) {
					sampleState = dynamics.getValue(node.getID());
				}
	
				else {
					sampleState = Integer.parseInt(splits[splits.length - 1]); // samples states (or priors) should
																				// eventually be specified in the XML
				}
				
				if ( node.getMetaData("location")!=null) {
					buf.append("[&");
					buf.append(dynamics.typeTraitInput.getName() + "="
							+ dynamics.getStringStateValue((int) node.getMetaData("location")));
					buf.append(']');
				}
	
			}
	
			buf.append(":");
			if (substitutions) {
				appendDouble(buf, (lastHeight-node.getHeight()) * branchRateModel.getRateForBranch(node));
			} else {
				appendDouble(buf, lastHeight-node.getHeight());
			}
		}
		return buf.toString();
	}


	public Node getRoot(){
		return mappedTree.getRoot();
	}
	
	@Override
	public void close(PrintStream out) {
		mappedTree.close(out);
	}

	public int getNodeState(int nr) {
		return (int) mappedTree.getNode(nr).getMetaData("location");
	}

}
