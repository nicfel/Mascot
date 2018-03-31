/*
 * Euler2ndOrderCPU.h
 *
 *  Created on: 30/03/2018
 *      Author: remco
 */

#ifndef EULER2NDORDERCPU_H_
#define EULER2NDORDERCPU_H_

#include "Euler2ndOrder.h"

class Euler2ndOrderCPU : public Euler2ndOrder {
	double epsilon;
	double max_step;

	double * migration_rates; // flattened square matrix of migration rates
	int n, n2;// dimension of migration rate matrix and indicators matrix
	int * multiplicator;
	int * indicators;
	double * coalescent_rates;
	double probs;
	int lineages;
	int states;
	int dimension;
	double * sumStates;
	bool hasIndicators;
	bool hasMultiplicator;
	double * tCR;
	double * sumDotStates;

	int iterations;

	double * linProbs_tmpdt;
	double * linProbs_tmpddt;
	double * linProbs_tmpdddt;

public:
	Euler2ndOrderCPU();
	virtual ~Euler2ndOrderCPU();
	void setup(int maxSize);
	void init(double * migration_rates, int n, double * coalescent_rates, int lineages, int states, double epsilon, double max_step);
	void initWithIndicators(double * migration_rates, int * indicators, double * coalescent_rates, int lineages, int states, double epsilon, double max_step);
	void calculateValues(double duration, double * p, int length);

private:
	void calculateValues(double duration, double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length);
	double updateP (double duration, double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length);
	double maxAbs(double * pDotDotDot, int length);
	void normalise(const int i, double * p);
	void updateP2(const double timeStep, const double timeStepSquare, double * p, const int length, double * pDot, double * pDotDot);
	void computeDerivatives (double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length);
	void calcSumStates(double  * sumStates, const double * p);
	void computeSecondDerivate (double * p, double * pDot, double * pDotDot, int length);
	void approximateThirdDerivate (double * pDotDot, double * pDotDotDot, int length);
	void computeDerivativesWithMultiplicator(double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length);
	void computeSecondDerivateWithMultiplicator(double * p, double * pDot, double * pDotDot, int length);

};

#endif /* EULER2NDORDERCPU_H_ */
