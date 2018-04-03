/*
 * Euler2ndOrder.cpp
 *
 *  Created on: 30/03/2018
 *      Author: remco
 */

#include "Euler2ndOrder.h"

Euler2ndOrder::Euler2ndOrder() {
	// TODO Auto-generated constructor stub

}

Euler2ndOrder::~Euler2ndOrder() {
	// TODO Auto-generated destructor stub
}

void Euler2ndOrder::setup(int maxSize, int states, double epsilon, double max_step) {

}

void Euler2ndOrder::init(double * migration_rates, int ratesCount, double * coalescent_rates, int lineages) {}
void Euler2ndOrder::initWithIndicators(double * migration_rates, int * indicators, double * coalescent_rates, int lineages) {}
void Euler2ndOrder::calculateValues(double duration, double * p, int length) {}
void Euler2ndOrder::setUpDynamics(int count, double * migration_rates, double * coalescent_rates, double * next_rate_shift) {}
