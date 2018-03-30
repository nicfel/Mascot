/*
 * Euler2ndOrder.h
 *
 *  Created on: 30/03/2018
 *      Author: remco
 */

#ifndef EULER2NDORDER_H_
#define EULER2NDORDER_H_

class Euler2ndOrder {
public:
	Euler2ndOrder();
	virtual ~Euler2ndOrder();
	void setup(int maxSize);
	void init(double * migration_rates, int n, double * coalescent_rates, int lineages, int states, double epsilon, double max_step);
	void initWithIndicators(double * migration_rates, int * indicators, double * coalescent_rates, int lineages, int states, double epsilon, double max_step);
	void calculateValues(double duration, double * p, int length);

};

#endif /* EULER2NDORDER_H_ */
