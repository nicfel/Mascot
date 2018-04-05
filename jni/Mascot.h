/*
 * Mascot.h
 *
 *  Created on: 5/04/2018
 *      Author: remco
 */

#ifndef MASCOT_H_
#define MASCOT_H_

#include "Euler2ndOrderCPU.h"

class Mascot {
public:
	 int states;
private:

	 int nrSamples;
	 double* stateProbabilities;

     int nrLineages;

    // current rates
     double* coalescentRates;


    // Set up for lineage state probabilities
     int * activeLineages;
     int activeLineagesCount;

	 double* linProbs;
	 double* linProbsNew;
	 int linProbsLength;

    // store the linProbs, multiplicators and logP's at coalescent points in jagged arrays from last time
     double* coalLinProbs;
     int * coalLinProbsLengths;
     double* coalLogP;
     int* coalRatesInterval;

    // deep store the things above for MCMC
     double* storeLinProbs;
     int * storedCoalLinProbsLengths;
     double* storeLogP;
     int* storeRatesInterval;

	 double * nextTreeEvents;
	 double * storedNextTreeEvents;
	 double * nextRateShifts;
	 double * storedNextRateShifts;

    // check if this is the first calculation
     bool first;


	// maximum integration error tolerance
     Euler2ndOrderCPU * euler;

     int * nodeType;

     int * parents;
     int * lineagesAdded;
     int * lineagesRemoved;
     double * intervals;

     int * storedLineagesAdded;
     int intervalCount, nodeCount;

     double * linProbs_tmp;

	 double * coalescentRatesCache;
	 double * migrationRatesCache;
	 int ** indicatorsCache;
	 double * nextRateShiftCache;
	 int rateShiftCount;

	 double * currentCoalescentRates;

	 double logP;
public:
	Mascot(int * nodeType, int states, double epsilon, double max_step, int sampleCount, int nodeCount, int intervalCount);
	virtual ~Mascot();
    double calculateLogP(bool dynamicsIsDirty, int firstDirtyInterval, int* lineagesAdded, int* lineagesRemoved, double* intervals, int* parents);
    void addActiveLineages(const int pos1, const int coalLines1);
    void addActiveLineages(int newLineage);
    void removeActiveLineageAt(int pos0);
    int indexOf(int value);


	double* getCoalescentRate(int i) ;
	double getRateShiftInterval(int i) ;
	void setUpDynamics(int count, double * migration_rates, double * coalescent_rates, double * next_rate_shift);
	double doEuler(double nextEventTime, int ratesInterval) ;
    void sample(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) ;
    double coalesce(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) ;
    void getStateProb(int nr, double * p) ;
    void getRootState(double * p) ;
    void storeNode(int storingTreeInterval, int storingRatesInterval, double* storeLinProbs,double probability, double nextTreeEvent, double nextRateShift);
    int restoreNode(int restoringInterval);
	void store() ;
    void storeLinP() ;
	void restore();
};

#endif /* MASCOT_H_ */
