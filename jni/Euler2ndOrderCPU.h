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
public:
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

	int rateShiftCount;
	double * migrationRatesCache;
	double * coalescentRatesCache;
	int * indicatorsRatesCache;
	//double * nextRateShiftCache;
public:
	Euler2ndOrderCPU();
	virtual ~Euler2ndOrderCPU();
	virtual void setup(int maxSize, int states, double epsilon, double max_step);
	virtual void init(double * migration_rates, int n, double * coalescent_rates, int lineages);
	virtual void initWithIndicators(double * migration_rates, int * indicators, double * coalescent_rates, int lineages);
	virtual void calculateValues(double duration, double * p, int length);
	virtual void setUpDynamics(int count, double * migration_rates, double * coalescent_rates, double * next_rate_shift);

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

#define sub(X) \
class X : public Euler2ndOrderCPU {\
public: X();\
virtual ~X();\
void setup(int maxSize, int states, double epsilon, double max_step);\
void init(double * migration_rates, int n, double * coalescent_rates, int lineages);\
void initWithIndicators(double * migration_rates, int * indicators, double * coalescent_rates, int lineages);\
void calculateValues(double duration, double * p, int length);\
void setUpDynamics(int count, double * migration_rates, double * coalescent_rates, double * next_rate_shift);\
\
private:\
void calculateValues(double duration, double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length);\
double updateP (double duration, double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length);\
double maxAbs(double * pDotDotDot, int length);\
void normalise(const int i, double * p);\
void updateP2(const double timeStep, const double timeStepSquare, double * p, const int length, double * pDot, double * pDotDot);\
void computeDerivatives (double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length);\
void calcSumStates(double  * sumStates, const double * p);\
void computeSecondDerivate (double * p, double * pDot, double * pDotDot, int length);\
void approximateThirdDerivate (double * pDotDot, double * pDotDotDot, int length);\
void computeDerivativesWithMultiplicator(double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length);\
void computeSecondDerivateWithMultiplicator(double * p, double * pDot, double * pDotDot, int length);\
\
};

sub(Euler2ndOrderCPU2);
sub(Euler2ndOrderCPU3);
sub(Euler2ndOrderCPU4);
sub(Euler2ndOrderCPU5);
sub(Euler2ndOrderCPU6);
sub(Euler2ndOrderCPU7);
sub(Euler2ndOrderCPU8);
sub(Euler2ndOrderCPU9);
sub(Euler2ndOrderCPU10);
sub(Euler2ndOrderCPU11);
sub(Euler2ndOrderCPU12);
sub(Euler2ndOrderCPU13);
sub(Euler2ndOrderCPU14);
sub(Euler2ndOrderCPU15);

#endif /* EULER2NDORDERCPU_H_ */
