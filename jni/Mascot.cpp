#include "Mascot.h"
#include <cmath>
#include <cfloat>
#include <iostream>
#include <immintrin.h>

#include <stdio.h>


inline double min(const double x, const double y) {return x < y ? x : y;}
inline double max(const double x, const double y) {return x > y ? x : y;}

inline void printArray(int * array, int n) {
//	return;
	fprintf(stderr, "[");
	for (int i = 0; i < n; i++) {
		fprintf(stderr, "%d", array[i]);
		if (i < n-1) {
			fprintf(stderr, " ");
		}
	}
	fprintf(stderr, "]\n");
}
inline void printArrayF(double * array, int n) {
//	return;
	fprintf(stderr, "[");
	for (int i = 0; i < n; i++) {
		fprintf(stderr, "%.2f", array[i]);
		if (i < n-1) {
			fprintf(stderr, " ");
		}
	}
	fprintf(stderr, "]\n");
}


inline void SystemArraycopy(double * src, const int start, double * dest, const int offset, const int count) {
	memcpy(dest + offset, src + start, count * sizeof(double));
}

inline void SystemArraycopy(double * src, double * dest, const int count) {
	memcpy(dest, src, count * sizeof(double));
}

inline void SystemArraycopyI(int * src, const int start, int * dest, const int offset, const int count) {
	memcpy(dest + offset, src + start, count * sizeof(int));
}

inline void SystemArraycopy(int * src, int * dest, const int count) {
	memcpy(dest, src, count * sizeof(int));
}

Mascot::Mascot(int * nodeType, int states, double epsilon, double max_step, int sampleCount, int nodeCount, int intervalCount, bool useCache) {
	(*this).nodeCount = nodeCount;
    	(*this).states = states;
    	(*this).nodeType = nodeType;
    	(*this).useCache = useCache;

    	fprintf(stderr, "nodeCount=%d states=%d sampleCount=%d intervalCount=%d useCache=%d\n", nodeCount, states, sampleCount, intervalCount, useCache);
    nrSamples = sampleCount + 1;
    	stateProbabilities = new double[nrSamples * states];


    	// initialize storing arrays and ArrayLists
    	storedLineagesAdded = new int[intervalCount];
    	coalLinProbs = new double[intervalCount * intervalCount * states];
    	storeLinProbs = new double[intervalCount * intervalCount * states];
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

    	activeLineages = new int[nodeCount + 1];
    	activeLineagesCount = 0;

    	int MAX_SIZE = intervalCount * states;
    	linProbs_tmp = new double[MAX_SIZE];
    	linProbs = new double[MAX_SIZE];
    	linProbsNew = new double[MAX_SIZE];

    	switch (states) {
    	case 2: euler = new Euler2ndOrderCPU2();break;
    	case 3: euler = new Euler2ndOrderCPU3();break;
    	case 4: euler = new Euler2ndOrderCPU4();break;
    	case 5: euler = new Euler2ndOrderCPU5();break;
    	case 6: euler = new Euler2ndOrderCPU6();break;
    	case 7: euler = new Euler2ndOrderCPU7();break;
    	case 8: euler = new Euler2ndOrderCPU8();break;
    	case 9: euler = new Euler2ndOrderCPU9();break;
    	case 10: euler = new Euler2ndOrderCPU10();break;
    	default: euler = new Euler2ndOrderCPU();
    	}
    	euler->setup(MAX_SIZE, states, epsilon, max_step);

    	// reserve memory for dynamics cache once we know how many rateshifts there are
    	// for now, set to 0
    	nextRateShiftCache = 0;
    	migrationRatesCache = 0;
    	coalescentRatesCache = 0;
    	//currentCoalescentRates = 0;
    	indicatorsCache = 0;

    	coalescentRates = 0; // == current Coalescent Rates

    	nrLineages = 0;
    	rateShiftCount = 0;
    	intervals = new double[intervalCount];
    	lineagesAdded = new int[intervalCount];
    	parents = new int[intervalCount];
    	lineagesRemoved = new int[intervalCount*2];
    	linProbsLength = 0;
    	(*this).intervalCount = intervalCount;
    	logP = 0;

    	first = true;
    	debug = false;
    	callCount = 0;

    }



double Mascot::calculateLogP(bool dynamicsIsDirty, int firstDirtyInterval, int* lineagesAddedIn, int* lineagesRemovedIn, double* intervalsIn, int* parentsIn) {
	SystemArraycopyI(lineagesAddedIn, 0, lineagesAdded, 0, intervalCount);
	SystemArraycopyI(lineagesRemovedIn, 0, lineagesRemoved, 0, intervalCount * 2);
	SystemArraycopy(intervalsIn, 0, intervals, 0, intervalCount);
	SystemArraycopyI(parentsIn, 0, parents, 0, nodeCount);

	//fprintf(stderr,"callCount = %d\n", callCount);
	//if (callCount == 244) {
	//	debug = true;
	//}
	//fprintf(stderr, "intervals"); printArrayF(intervals, intervalCount);

	//printArray(lineagesRemoved, 20);

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

      if (useCache && !first && !dynamicsIsDirty && firstDirtyInterval > 2) {
        // restore the likelihood to last known good place
    	  int pos0 = -1, pos1 = -1;
    	  do {

        	nextEventTime = nextTreeEvent;

    		// Check if the last interval was reached
    		if (treeInterval == intervalCount){
    			logP = coalLogP[intervalCount - 1];
    			return logP;
    		}

    		bool isDirty;
    		if (lineagesRemoved[treeInterval*2] >= 0) { // == IntervalType.COALESCENT) {
 	           	int coalLines0 = lineagesRemoved[treeInterval * 2 + 0];
 	           	int coalLines1 = lineagesRemoved[treeInterval * 2 + 1];
 	           	pos0 = indexOf(coalLines0);
 	           	pos1 = indexOf(coalLines1);
 	           	if (pos0 < 0 || pos1 < 0) {
 	           		std::cerr << coalLines0 <<  " " << coalLines1 <<  " " <<activeLineages << std::endl;
 	           		std::cerr << "daughter lineages at coalescent event not found" << std::endl;
 	           		std::cerr << "coalesceX went wrong at 1"<< std::endl;
 	           		exit(1);
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

        } while (nextTreeEvent <= INFINITY);
    }


		coalescentRates = getCoalescentRate(ratesInterval);
		nrLineages = activeLineagesCount;
		linProbsLength = nrLineages * states;

        // Calculate the likelihood
        do {
        	nextEventTime = min(nextTreeEvent, nextRateShift);
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
//std::cerr << "logP = " << logP << " " << treeInterval << " " << ratesInterval << std::endl;
        	if (logP == -INFINITY) {
        		return logP;
        	}
        } while (nextTreeEvent <= INFINITY);
//std::cerr << "exiting with logP = " << logP << std::endl;

        first = false;
        callCount++;
		return logP;
    }


	void Mascot::addActiveLineages(const int pos1, const int coalLines1) {
		SystemArraycopyI(activeLineages, pos1, activeLineages, pos1 + 1, activeLineagesCount - pos1);
		activeLineages[pos1] = coalLines1;
		activeLineagesCount++;
		//printArray(activeLineages, activeLineagesCount);
	}

	void Mascot::addActiveLineages(int newLineage) {
		activeLineages[activeLineagesCount++] = newLineage;
		//printArray(activeLineages, activeLineagesCount);
	}

	void  Mascot::removeActiveLineageAt(int pos0) {
		SystemArraycopyI(activeLineages, pos0 + 1, activeLineages, pos0, activeLineagesCount - pos0);
		activeLineagesCount--;
		//printArray(activeLineages, activeLineagesCount);
	}

	int  Mascot::indexOf(int value) {
		for (int i = 0; i < activeLineagesCount; i++) {
			if (activeLineages[i] == value) {
				return i;
			}
		}
		return -1;
	}

	double* Mascot::getCoalescentRate(int i) {
		if (i >= rateShiftCount) {
			return coalescentRatesCache + (rateShiftCount-1) * states;
		} else {
			return coalescentRatesCache + i * states;
		}
	}

	double* Mascot::getMigrationRates(int i) {
		if (i >= rateShiftCount) {
			return migrationRatesCache +  (rateShiftCount-1) * states * states;
		} else {
			return migrationRatesCache +  i * states * states;
		}
	}

	double Mascot::getRateShiftInterval(int i) {
    		if (i >= rateShiftCount) {
    			return INFINITY;
    		} else {
			return nextRateShiftCache[i];
    		}
	}

	void Mascot::setUpDynamics(int count, double * coalescent_rates, double * migration_rates, double * next_rate_shift) {
		if (rateShiftCount != count) {
			rateShiftCount = count;
			coalescentRatesCache = new double[rateShiftCount * states];
			migrationRatesCache = new double[rateShiftCount * states * states];
        		nextRateShiftCache = new double[rateShiftCount];
		}
		SystemArraycopy(coalescent_rates, 0, coalescentRatesCache, 0, rateShiftCount * states);
		SystemArraycopy(migration_rates, 0, migrationRatesCache, 0, rateShiftCount * states * states);
		SystemArraycopy(next_rate_shift, 0, nextRateShiftCache, 0, rateShiftCount);
	}


	double Mascot::doEuler(double nextEventTime, int ratesInterval) {
		SystemArraycopy(linProbs,0,linProbs_tmp,0,linProbsLength);
		linProbs_tmp[linProbsLength] = 0;
		linProbs[linProbsLength-1] = 0;

		if (debug) {
			fprintf(stderr,"%f\n", nextEventTime);
			fprintf(stderr,"caol ");
			printArrayF(coalescentRates, states);
			fprintf(stderr,"imgr ");
			printArrayF(getMigrationRates(ratesInterval), states * states);
			fprintf(stderr,"p ");
			printArrayF(linProbs_tmp, linProbsLength);
		}
		euler->init(getMigrationRates(ratesInterval),
				states * states,
				coalescentRates,
				nrLineages);
		euler->calculateValues(nextEventTime, linProbs_tmp, linProbsLength + 1);

		SystemArraycopy(linProbs_tmp,0,linProbs,0,linProbsLength);

		return linProbs_tmp[linProbsLength];
	}



    void Mascot::sample(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
			if (debug) {
				fprintf(stderr,"sample activeLineages %d = ", currTreeInterval);
				printArray(activeLineages, activeLineagesCount);
			}
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

    double Mascot::coalesce(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
    	if (debug) {
    		fprintf(stderr,"coalesce activeLineages %d = ", currTreeInterval);
    		printArray(activeLineages, activeLineagesCount);
    	}

    	int coalLines0 = lineagesRemoved[currTreeInterval * 2 + 0];
    	int coalLines1 = lineagesRemoved[currTreeInterval * 2 + 1];

    	const int daughterIndex1 = indexOf(coalLines0);//.getNr());
		const int daughterIndex2 = indexOf(coalLines1);//.getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			std::cerr << coalLines0/*.getNr()*/ << " " << coalLines1/*.getNr()*/ << " " << activeLineages << std::endl;
			std::cerr << "daughter lineages at coalescent event not found" << std::endl;
			exit(1);
		}

		int offset = (parents[coalLines0] - nrSamples) * states;
		//double * lambda = stateProbabilities[parents[coalLines0] - nrSamples];

//		std::cerr << coalLines0/*.getNr()*/ << " " << coalLines1/*.getNr()*/ << " " << activeLineages  << " " << offset << std::endl;
		/*
		 * Calculate the overall probability for two strains to coalesce
		 * independent of the state at which this coalescent event is
		 * supposed to happen
		 */
        for (int k = 0; k < states; k++) {
        		double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];

        		//std::cerr << pairCoalRate << "=" << coalescentRates[k] << "*" << linProbs[daughterIndex1*states + k] << "*" << linProbs[daughterIndex2*states + k] << std::endl;
			if (!(pairCoalRate == NAN)){
				stateProbabilities[offset + k] = pairCoalRate;
			} else {
				return -INFINITY;
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
//		std::cerr << "Sum is: " << sum << std::endl;

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


		double min = INFINITY;
        for (int i = 0; i < states; i++) {
			if (stateProbabilities[offset + i] < min) {min = stateProbabilities[offset + i];}
		}
		if (min < 0.0){
			std::cerr << "Coalescent probability is: " << min << std::endl;
			return -INFINITY;
		}

		if (sum==0) {
			return -INFINITY;
		}

		// store the node
		sum = log(sum);
		storeNode(currTreeInterval, currRatesInterval, linProbs, logP + sum, nextTreeEvent, nextRateShift);

//		std::cerr << "Final coalescent probability is: " << sum << std::endl;

//		if (sum==0)
//			return -INFINITY;
//		else
		return sum;
    }


    void Mascot::getStateProb(int nr, double * p) {
    	for (int i = 0; i < states; i++) {
    		p[i] = stateProbabilities[(nr - nrSamples) * states + i];
    	}
    }

    void Mascot::getRootState(double * p) {
    	for (int i = 0; i < states; i++) {
    		p[i] = stateProbabilities[(nrSamples-2) * states + i];
    	}
    }

    void Mascot::storeNode(int storingTreeInterval, int storingRatesInterval, double* storeLinProbs,
    		double probability, double nextTreeEvent, double nextRateShift) {
    	if (!useCache) {
    		return;
    	}
    	coalRatesInterval[storingTreeInterval] = storingRatesInterval;
    	int offset = 0;
    	if (storingTreeInterval > 0) {
    		offset = coalLinProbsLengths[storingTreeInterval-1];
    	}
    	SystemArraycopy(storeLinProbs, 0, coalLinProbs, offset, linProbsLength);
    	coalLinProbsLengths[storingTreeInterval] = offset + linProbsLength;
    	coalLogP[storingTreeInterval] = probability;
    	nextTreeEvents[storingTreeInterval] = nextTreeEvent;
    	nextRateShifts[storingTreeInterval] = nextRateShift;
    }

    int Mascot::restoreNode(int restoringInterval){
    	//Log.warning("Restoring " + first + " " + restoringInterval);
    	int offset = 0;
    	if (restoringInterval > 0) {
    		offset = coalLinProbsLengths[restoringInterval-1];
    	}
    	linProbsLength = coalLinProbsLengths[restoringInterval] - offset;
    	SystemArraycopy(coalLinProbs, offset, linProbs, 0, linProbsLength);

    	logP = coalLogP[restoringInterval];
    	return coalRatesInterval[restoringInterval + 1];

    }

	void Mascot::store() {
    	if (!useCache) {
    		return;
    	}
    	storeLinP();
    	SystemArraycopy(coalLogP, 0, storeLogP, 0, intervalCount);
    	SystemArraycopyI(coalRatesInterval, 0, storeRatesInterval, 0, intervalCount);


    	SystemArraycopy(nextTreeEvents, 0, storedNextTreeEvents, 0, intervalCount);
    	SystemArraycopy(nextRateShifts, 0, storedNextRateShifts, 0, intervalCount);

    	SystemArraycopyI(lineagesAdded, 0, storedLineagesAdded, 0, intervalCount);

    }

    void Mascot::storeLinP() {
    	SystemArraycopyI(coalLinProbsLengths, 0, storedCoalLinProbsLengths, 0, intervalCount);
    	SystemArraycopy(coalLinProbs, 0, storeLinProbs, 0, coalLinProbsLengths[intervalCount - 1]);
	}


	void Mascot::restore(){
    	if (!useCache) {
    		return;
    	}
    	// restore intermediate results
    	double * tmp = storeLogP;
    	storeLogP = coalLogP;
    	coalLogP = tmp;

    	tmp = coalLinProbs;
    	coalLinProbs = storeLinProbs;
    	storeLinProbs = tmp;

    	int * tmp2 = coalLinProbsLengths;
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

		tmp2 = lineagesAdded;
		lineagesAdded = storedLineagesAdded;
		storedLineagesAdded = tmp2;

    }



	Mascot::~Mascot() {
		// TODO Auto-generated destructor stub
	}

