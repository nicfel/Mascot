/*
 * Euler2ndOrder.h
 *
 *  Created on: 30/03/2018
 *      Author: remco
 */

#ifndef EULER2NDORDER_H_
#define EULER2NDORDER_H_

#include <string.h>
#include <stdio.h>

class Euler2ndOrder {
public:
	Euler2ndOrder();
	virtual ~Euler2ndOrder();
	void setup(int maxSize, int states, double epsilon, double max_step);
	void init(double * migration_rates, int n, double * coalescent_rates, int lineages);
	void initWithIndicators(double * migration_rates, int * indicators, double * coalescent_rates, int lineages);
	void calculateValues(double duration, double * p, int length);
	void setUpDynamics(int count, double * migration_rates, double * coalescent_rates, double * next_rate_shift);
};

#endif /* EULER2NDORDER_H_ */
