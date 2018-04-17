#include "Euler2ndOrderCPU.h"
#include <cmath>
#include <cfloat>
#include <iostream>
#include <immintrin.h>
using namespace std;


inline double min(const double x, const double y) {return x < y ? x : y;}
inline double max(const double x, const double y) {return x > y ? x : y;}

Euler2ndOrderCPU8::Euler2ndOrderCPU8() {
	fprintf(stderr, "Creating Euler2ndOrderCPU8\n");
	max_step = 0 , epsilon = 0;
	migration_rates = NULL; // flattened square matrix of migration rates
	n = 0, n2 = 0;// dimension of migration rate matrix and indicators matrix
	multiplicator = NULL;
	indicators = NULL;
	coalescent_rates = NULL;
	probs = 0;
	lineages = 0;
	states = 0;
	dimension = 0;
	sumStates  = NULL;
	hasIndicators = false;
	hasMultiplicator = false;
	tCR  = NULL;
	sumDotStates = NULL;

	iterations = 0;

	linProbs_tmpdt = NULL;
	linProbs_tmpddt = NULL;
	linProbs_tmpdddt = NULL;

	rateShiftCount = 0;

	migrationRatesCache = NULL;
	coalescentRatesCache = NULL;
	indicatorsRatesCache = NULL;
	//nextRateShiftCache = NULL;
}

Euler2ndOrderCPU8::~Euler2ndOrderCPU8() {
}



void Euler2ndOrderCPU8::init(double * migration_rates, int rateCount, double * coalescent_rates, int lineages) {

//	fprintf(stderr,"x");
		(*this).migration_rates = migration_rates;
		(*this).n = (int)(sqrt(rateCount) + 0.5);
		(*this).coalescent_rates = coalescent_rates;
		(*this).lineages = lineages;
		(*this).dimension = (*this).lineages*(*this).states;
		hasIndicators = false;
		hasMultiplicator = false;

		iterations=0;


//		int k = 0;
//		for (int i = 0;  i < n; i++) {
//			fprintf(stderr,"\n");
//			for (int j = 0; j < n; j++) {
//				fprintf(stderr,"%f ", migration_rates[k++]);
//			}
//			fprintf(stderr,"\n");
//		}
	}

void Euler2ndOrderCPU8::initWithIndicators(double * migration_rates, int * indicators, double * coalescent_rates, int lineages) {
//		(*this).max_step = max_step;
//		(*this).epsilon = epsilon;
//		(*this).migration_rates = migration_rates;
//		n = (int)(sqrt(migration_rates.length) + 0.5);
//		(*this).indicators = indicators;
//		n2 = indicators.length / 2;
//		(*this).coalescent_rates = coalescent_rates;
//		(*this).lineages = lineages;
//		(*this).states = states;
//		(*this).dimension = (*this).lineages*(*this).states;
//		sumStates = new double[states];tCR = new double[states];
//		hasIndicators = true;
//		hasMultiplicator = false;
//
//		iterations=0;
	}


void Euler2ndOrderCPU8::setup(int maxSize, int states, double epsilon, double max_step) {
		linProbs_tmpdt = new double[maxSize];
		linProbs_tmpddt = new double[maxSize];
		linProbs_tmpdddt = new double[maxSize];
		memset(linProbs_tmpdt, 0, maxSize * sizeof(double));
		memset(linProbs_tmpddt, 0, maxSize * sizeof(double));
		memset(linProbs_tmpdddt, 0, maxSize * sizeof(double));
		fprintf(stderr,"Reserving %d\n", maxSize);
		(*this).states = states;
		(*this).epsilon = epsilon;
		(*this).max_step = max_step;

		sumStates = new double[states];
		tCR = new double[states];
		sumDotStates = new double[states];
}


void Euler2ndOrderCPU8::setUpDynamics(int count, double * migration_rates, double * coalescent_rates, double * next_rate_shift) {
	if (	rateShiftCount == 0) {
		(*this).rateShiftCount = count;
		migrationRatesCache = new double[rateShiftCount * states * states];
		coalescentRatesCache = new double[rateShiftCount * states];
		indicatorsRatesCache = new int[rateShiftCount * states];
		//nextRateShiftCache = new double[rateShiftCount];
	}
	memcpy(migrationRatesCache, migration_rates, rateShiftCount * states * states * sizeof(double));
	memcpy(coalescentRatesCache, coalescent_rates, rateShiftCount * states * sizeof(double));
	//memcpy(nextRateShiftCache, next_rate_shift, rateShiftCount * sizeof(double));
}

void Euler2ndOrderCPU8::calculateValues(double duration, double * p, int length) {
		double * pDot = linProbs_tmpdt;
		double * pDotDot = linProbs_tmpddt;
		double * pDotDotDot = linProbs_tmpdddt;
		calculateValues(duration, p, pDot, pDotDot, pDotDotDot, length);
}

//void Euler2ndOrderCPU8::calculateValues(double duration, double * p, double * pDot, double * pDotDot, double * pDotDotDot) {
//		calculateValues(duration, p, pDot, pDotDot, pDotDotDot, pDot.length);
//	}

void Euler2ndOrderCPU8::calculateValues(double duration, double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length) {
//	fprintf(stderr,"states=%d eps=%f max_step=%f\n", states, epsilon, max_step);
//	printArray(coalescent_rates, states);
//	printArray(migration_rates, states);
		memset(pDotDot, 0, length * sizeof(double));
		memset(pDotDotDot, 0, length * sizeof(double));

		if (hasMultiplicator) {
			while (duration > 0) {
				iterations++;
				//pDot = new double[length];
				memset(pDot, 0, length * sizeof(double));
				computeDerivativesWithMultiplicator(p, pDot, pDotDot, pDotDotDot, length);
				computeSecondDerivateWithMultiplicator(p, pDot, pDotDot, length);
				approximateThirdDerivate(pDotDot, pDotDotDot, length);
				duration = updateP(duration, p, pDot, pDotDot, pDotDotDot, length - 1);

				if (iterations>10000) {
					std::cerr << "too many iterations, return negative infinity" << std::endl;
					p[length-1] = DBL_MIN;//Double.NEGATIVE_INFINITY;
					break;
				}
			}
		} else {
			while (duration > 0) {
				iterations++;
				//pDot = new double[length];
				memset(pDot, 0, length * sizeof(double));
				computeDerivatives(p, pDot, pDotDot, pDotDotDot, length);
				computeSecondDerivate(p, pDot, pDotDot, length);
				approximateThirdDerivate(pDotDot, pDotDotDot, length);
				duration = updateP(duration, p, pDot, pDotDot, pDotDotDot, length - 1);

				if (iterations>10000) {
					std::cerr << "too many iterations, return negative infinity" << std::endl;
					p[length-1] = DBL_MIN;
					break;
				}
			}
		}
	}



double Euler2ndOrderCPU8::updateP (double duration, double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length) {
		const double max_dotdotdot = maxAbs(pDotDotDot, length);

		//double timeStep = FastMath.min(FastMath.pow(epsilon*6/max_dotdotdot, C), FastMath.min(duration, max_step));

		double timeStep = min(cbrt(epsilon*6/max_dotdotdot), min(duration, max_step));
		double timeStepSquare = timeStep * timeStep * 0.5;

		for (int i = 0; i < length; i++) {
			double new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
			double diff = abs(new_val - p[i]);
			while (new_val > 1 || new_val < 0 || diff > 0.2) {
				timeStep *= 0.9;
				timeStepSquare = timeStep * timeStep * 0.5;
				new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
				diff = abs(new_val - p[i]);
			}
		}

		updateP2(timeStep, timeStepSquare, p, length + 1, pDot, pDotDot);

		// normalize to ensure stability
		for (int i = 0; i < lineages; i ++) {
			normalise(i, p);
		}

		duration -= timeStep;
		return duration;
	}



double Euler2ndOrderCPU8::maxAbs(double * pDotDotDot, int length) {
		double max_dotdotdot = 0.0;
		for (int i = 0; i < length; i++) {
			max_dotdotdot = max(max_dotdotdot, abs(pDotDotDot[i]));
		}
		return max_dotdotdot;
	}

void Euler2ndOrderCPU8::normalise(const int i, double * p) {
		const int k = states * i;
		double *p2 = p + k;

		double linSum = 0;
		
{
			linSum += *p2++;
		}
{
			linSum += *p2++;
		}
{
			linSum += *p2++;
		}
{
			linSum += *p2++;
		}
{
			linSum += *p2++;
		}
{
			linSum += *p2++;
		}
{
			linSum += *p2++;
		}
{
			linSum += *p2++;
		}
		p2 = p + k;
{
			*p2++ /= linSum;
		}
{
			*p2++ /= linSum;
		}
{
			*p2++ /= linSum;
		}
{
			*p2++ /= linSum;
		}
{
			*p2++ /= linSum;
		}
{
			*p2++ /= linSum;
		}
{
			*p2++ /= linSum;
		}
{
			*p2++ /= linSum;
		}
	}

void Euler2ndOrderCPU8::updateP2(const double timeStep, const double timeStepSquare, double * p, const int length, double * pDot,
			double * pDotDot) {
		for (int i = 0; i < length; i++) {
			p[i] += pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
		}
	}

void Euler2ndOrderCPU8::computeDerivatives (double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length) {

		double migrates;
		// Compute the sum of line state probabilities for each state
		memset(sumStates, 0, states * sizeof(double));
		calcSumStates(sumStates, p);

		// Calculate the change in the lineage state probabilities for every lineage in every state
		int currlin = 0, j;
		for (int i = 0; i<lineages; i++) {

			double sumCoal = 0;
			int k = currlin;
j = 0;
{
				tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
				sumCoal += p[k] * tCR[j];
				k++;
			}
j++;
{
				tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
				sumCoal += p[k] * tCR[j];
				k++;
			}
j++;
{
				tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
				sumCoal += p[k] * tCR[j];
				k++;
			}
j++;
{
				tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
				sumCoal += p[k] * tCR[j];
				k++;
			}
j++;
{
				tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
				sumCoal += p[k] * tCR[j];
				k++;
			}
j++;
{
				tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
				sumCoal += p[k] * tCR[j];
				k++;
			}
j++;
{
				tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
				sumCoal += p[k] * tCR[j];
				k++;
			}
j++;
{
				tCR[j] = coalescent_rates[j] * (sumStates[j] - p[k]);
				sumCoal += p[k] * tCR[j];
				k++;
			}
			pDot[length-1] -= sumCoal;

			k = currlin;
j = 0;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[k] = coal;
				pDotDotDot[k] = coal;
				pDot[k] += p[k] * coal;
				k++;
			} // j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[k] = coal;
				pDotDotDot[k] = coal;
				pDot[k] += p[k] * coal;
				k++;
			} // j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[k] = coal;
				pDotDotDot[k] = coal;
				pDot[k] += p[k] * coal;
				k++;
			} // j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[k] = coal;
				pDotDotDot[k] = coal;
				pDot[k] += p[k] * coal;
				k++;
			} // j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[k] = coal;
				pDotDotDot[k] = coal;
				pDot[k] += p[k] * coal;
				k++;
			} // j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[k] = coal;
				pDotDotDot[k] = coal;
				pDot[k] += p[k] * coal;
				k++;
			} // j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[k] = coal;
				pDotDotDot[k] = coal;
				pDot[k] += p[k] * coal;
				k++;
			} // j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[k] = coal;
				pDotDotDot[k] = coal;
				pDot[k] += p[k] * coal;
				k++;
			} // j
			currlin += states;
		}

		// Calculate the probability of a lineage changing states
//		if (hasIndicators) {
//			for (int j = 0; j < indicators.length/2; j++) {
//				int source = indicators[j * n2 + 0];
//				int sink = indicators[j * n2 + 1];
//				double mrate = migration_rates[source * n + sink];
//				int k = source;
//				int m = sink;
//				for (int i = 0; i<lineages; i++) {
//					migrates = p[k] * mrate;
//					pDot[m] += migrates;
//					pDot[k] -= migrates;
//					k += states;
//					m += states;
//				}
//			}
//		} else {
			int u = 0;
			for (int i = 0; i < lineages; i++) {
				// Calculate the probability of a lineage changing states
j = 0;
{
					int v = u;
					double pj = p[u];
					for (int k = j + 1; k < states; k++) {
						v++;
						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[v] * migration_rates[k * n + j] -
						pj * migration_rates[j * n + k];
						pDot[u] += migrates;
						pDot[v] -= migrates;
					} // j XXX
					u++;
				} // j
j++;
{
					int v = u;
					double pj = p[u];
					for (int k = j + 1; k < states; k++) {
						v++;
						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[v] * migration_rates[k * n + j] -
						pj * migration_rates[j * n + k];
						pDot[u] += migrates;
						pDot[v] -= migrates;
					} // j XXX
					u++;
				} // j
j++;
{
					int v = u;
					double pj = p[u];
					for (int k = j + 1; k < states; k++) {
						v++;
						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[v] * migration_rates[k * n + j] -
						pj * migration_rates[j * n + k];
						pDot[u] += migrates;
						pDot[v] -= migrates;
					} // j XXX
					u++;
				} // j
j++;
{
					int v = u;
					double pj = p[u];
					for (int k = j + 1; k < states; k++) {
						v++;
						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[v] * migration_rates[k * n + j] -
						pj * migration_rates[j * n + k];
						pDot[u] += migrates;
						pDot[v] -= migrates;
					} // j XXX
					u++;
				} // j
j++;
{
					int v = u;
					double pj = p[u];
					for (int k = j + 1; k < states; k++) {
						v++;
						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[v] * migration_rates[k * n + j] -
						pj * migration_rates[j * n + k];
						pDot[u] += migrates;
						pDot[v] -= migrates;
					} // j XXX
					u++;
				} // j
j++;
{
					int v = u;
					double pj = p[u];
					for (int k = j + 1; k < states; k++) {
						v++;
						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[v] * migration_rates[k * n + j] -
						pj * migration_rates[j * n + k];
						pDot[u] += migrates;
						pDot[v] -= migrates;
					} // j XXX
					u++;
				} // j
j++;
{
					int v = u;
					double pj = p[u];
					for (int k = j + 1; k < states; k++) {
						v++;
						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[v] * migration_rates[k * n + j] -
						pj * migration_rates[j * n + k];
						pDot[u] += migrates;
						pDot[v] -= migrates;
					} // j XXX
					u++;
				} // j
j++;
{
					int v = u;
					double pj = p[u];
					for (int k = j + 1; k < states; k++) {
						v++;
						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[v] * migration_rates[k * n + j] -
						pj * migration_rates[j * n + k];
						pDot[u] += migrates;
						pDot[v] -= migrates;
					} // j XXX
					u++;
				} // j
			} // lineages       		
//		}

		pDot[length-1] /= 2;

	}

void Euler2ndOrderCPU8::calcSumStates(double  * sumStates, const double * p) {
		
		double * ss;
		for (int i = 0; i < lineages; i++) {
			ss = sumStates;
{
				*ss++ += *p++;
			}
{
				*ss++ += *p++;
			}
{
				*ss++ += *p++;
			}
{
				*ss++ += *p++;
			}
{
				*ss++ += *p++;
			}
{
				*ss++ += *p++;
			}
{
				*ss++ += *p++;
			}
{
				*ss++ += *p++;
			}
		}
	}

void Euler2ndOrderCPU8::computeSecondDerivate (double * p, double * pDot, double * pDotDot, int length) {
		memset(sumDotStates, 0, states * sizeof(double));
		calcSumStates(sumDotStates, pDot);

		// Calculate the change in the lineage state probabilities for every lineage in every state
		int currlin = 0, j;
		for (int i = 0; i < lineages; i++) {
			double pCoalRate = 0.0;
			int k = currlin;
j = 0;
{
				pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
				k++;
			}
j++;
{
				pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
				k++;
			}
j++;
{
				pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
				k++;
			}
j++;
{
				pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
				k++;
			}
j++;
{
				pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
				k++;
			}
j++;
{
				pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
				k++;
			}
j++;
{
				pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
				k++;
			}
j++;
{
				pCoalRate += coalescent_rates[j] * (pDot[k] * (sumStates[j] - 2 * p[k]) + p[k] * (sumDotStates[j]));
				k++;
			}

			k = currlin;
j = 0;
{
				pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
				k++;
			}    			// j
j++;
{
				pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
				k++;
			}    			// j
j++;
{
				pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
				k++;
			}    			// j
j++;
{
				pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
				k++;
			}    			// j
j++;
{
				pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
				k++;
			}    			// j
j++;
{
				pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
				k++;
			}    			// j
j++;
{
				pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
				k++;
			}    			// j
j++;
{
				pDotDot[k] = pDotDot[k] * pDot[k] + p[k] * (pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[k]));
				k++;
			}    			// j

			pDotDot[length-1] -= pCoalRate;

			currlin += states;
		}    			// lineages 

		double migrates;

						// Calculate the probability of a lineage changing states
//		if (hasIndicators) {
//			for (int j = 0; j < indicators.length/2; j++) {
//				int source = indicators[j * n2 + 0];
//				int sink = indicators[j * n2 + 1];
//				double mrate = migration_rates[source * n + sink];
//				for (int i = 0; i<lineages; i++) {
//					migrates = pDot[states*i+source]*mrate;
//					pDotDot[states*i+sink] += migrates;
//					pDotDot[states*i+source] -= migrates;
//				}
//			}
//		} else {
			int u = 0;
			for (int i = 0; i<lineages; i++) {
				// Calculate the probability of a lineage changing states
j = 0;
{
					double pj = pDot[u];
					int v = u + 1;
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[v]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[u] += migrates;
						pDotDot[v] -= migrates;
						v++;
					}    			// j    	XXX
					u++;
				}    			// j
j++;
{
					double pj = pDot[u];
					int v = u + 1;
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[v]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[u] += migrates;
						pDotDot[v] -= migrates;
						v++;
					}    			// j    	XXX
					u++;
				}    			// j
j++;
{
					double pj = pDot[u];
					int v = u + 1;
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[v]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[u] += migrates;
						pDotDot[v] -= migrates;
						v++;
					}    			// j    	XXX
					u++;
				}    			// j
j++;
{
					double pj = pDot[u];
					int v = u + 1;
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[v]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[u] += migrates;
						pDotDot[v] -= migrates;
						v++;
					}    			// j    	XXX
					u++;
				}    			// j
j++;
{
					double pj = pDot[u];
					int v = u + 1;
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[v]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[u] += migrates;
						pDotDot[v] -= migrates;
						v++;
					}    			// j    	XXX
					u++;
				}    			// j
j++;
{
					double pj = pDot[u];
					int v = u + 1;
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[v]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[u] += migrates;
						pDotDot[v] -= migrates;
						v++;
					}    			// j    	XXX
					u++;
				}    			// j
j++;
{
					double pj = pDot[u];
					int v = u + 1;
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[v]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[u] += migrates;
						pDotDot[v] -= migrates;
						v++;
					}    			// j    	XXX
					u++;
				}    			// j
j++;
{
					double pj = pDot[u];
					int v = u + 1;
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[v]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[u] += migrates;
						pDotDot[v] -= migrates;
						v++;
					}    			// j    	XXX
					u++;
				}    			// j
			}    			// lineages    

//		}
		pDotDot[length-1] /= 2;
	}

void Euler2ndOrderCPU8::approximateThirdDerivate (double * pDotDot, double * pDotDotDot, int length) {
		double migrates;

		// Calculate the change in the lineage state probabilities for every lineage in every state
		for (int u = 0; u < length - 1; u++) {
			pDotDotDot[u] *= pDotDot[u];
		}

		// Calculate the probability of a lineage changing states
//		if (hasIndicators) {
//			for (int j = 0; j < indicators.length/2; j++) {
//				int source = indicators[j * n2 + 0];
//				int sink = indicators[j * n2 + 1];
//				double mrate = migration_rates[source * n + sink];
//				for (int i = 0; i<lineages; i++) {
//					migrates = pDotDot[states * i + source] * mrate;
//					pDotDotDot[states * i + sink] += migrates;
//					pDotDotDot[states * i + source] -= migrates;
//				}
//			}
//		} else {
		int k;
			for (int j = 0; j < states; j++) {
k = 0;
{
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
					for (int i = 0; i<lineages; i++) {
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
					} // XXX

				}
k++;
{
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
					for (int i = 0; i<lineages; i++) {
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
					} // XXX

				}
k++;
{
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
					for (int i = 0; i<lineages; i++) {
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
					} // XXX

				}
k++;
{
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
					for (int i = 0; i<lineages; i++) {
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
					} // XXX

				}
k++;
{
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
					for (int i = 0; i<lineages; i++) {
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
					} // XXX

				}
k++;
{
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
					for (int i = 0; i<lineages; i++) {
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
					} // XXX

				}
k++;
{
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
					for (int i = 0; i<lineages; i++) {
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
					} // XXX

				}
k++;
{
					double mrate = migration_rates[j * n + k];
					int u = j;
					int v = k;
					for (int i = 0; i<lineages; i++) {
						migrates = pDotDot[u] * mrate;
						pDotDotDot[v] += migrates;
						pDotDotDot[u] -= migrates;
						u += states;
						v += states;
					} // XXX

				}
			}
//		}
	}

void Euler2ndOrderCPU8::computeDerivativesWithMultiplicator(double * p, double * pDot, double * pDotDot, double * pDotDotDot, int length) {

		double migrates;
		int j;
		// Compute the sum of line state probabilities for each state
		memset(sumStates, 0, states * sizeof(double));
		{
			int u = 0;
			for (int i = 0; i<lineages; i++) {
				for (j = 0; j<states; j++) {
					sumStates[j] += multiplicator[i]*p[u++];
				}
			}
		}

		// Calculate the change in the lineage state probabilities for every lineage in every state
		//double * tCR =  new double[states];
		for (int i = 0; i<lineages; i++) {
			double sumCoal = 0;
			int currlin = states*i;
			for (j = 0; j<states; j++) {
				tCR[j] = coalescent_rates[j] * (sumStates[j] - p[currlin+j]);
				sumCoal += p[currlin+j]*tCR[j];
			}
			pDot[length-1] -= multiplicator[i]*sumCoal;
j = 0;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[currlin+j] = coal;
				pDotDotDot[currlin+j] = coal;
				pDot[currlin+j] += p[currlin+j] * coal;
			}		// j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[currlin+j] = coal;
				pDotDotDot[currlin+j] = coal;
				pDot[currlin+j] += p[currlin+j] * coal;
			}		// j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[currlin+j] = coal;
				pDotDotDot[currlin+j] = coal;
				pDot[currlin+j] += p[currlin+j] * coal;
			}		// j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[currlin+j] = coal;
				pDotDotDot[currlin+j] = coal;
				pDot[currlin+j] += p[currlin+j] * coal;
			}		// j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[currlin+j] = coal;
				pDotDotDot[currlin+j] = coal;
				pDot[currlin+j] += p[currlin+j] * coal;
			}		// j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[currlin+j] = coal;
				pDotDotDot[currlin+j] = coal;
				pDot[currlin+j] += p[currlin+j] * coal;
			}		// j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[currlin+j] = coal;
				pDotDotDot[currlin+j] = coal;
				pDot[currlin+j] += p[currlin+j] * coal;
			}		// j
j++;
{
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
				pDotDot[currlin+j] = coal;
				pDotDotDot[currlin+j] = coal;
				pDot[currlin+j] += p[currlin+j] * coal;
			}		// j

		}

		// Calculate the probability of a lineage changing states
//		if (hasIndicators) {
//			for (int j = 0; j < indicators.length/2; j++) {
//				int source = indicators[j * n2 + 0];
//				int sink = indicators[j * n2 + 1];
//				double mrate = migration_rates[source * n + sink];
//				for (int i = 0; i<lineages; i++) {
//					migrates = p[states*i+source]*mrate;
//					pDot[states*i+sink] += migrates;
//					pDot[states*i+source] -= migrates;
//				}
//			}
//		} else {
			for (int i = 0; i<lineages; i++) {
				int currlin = states*i;
				// Calculate the probability of a lineage changing states
j = 0;
{
					double pj = p[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDot[currlin+j] += migrates;
						pDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = p[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDot[currlin+j] += migrates;
						pDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = p[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDot[currlin+j] += migrates;
						pDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = p[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDot[currlin+j] += migrates;
						pDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = p[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDot[currlin+j] += migrates;
						pDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = p[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDot[currlin+j] += migrates;
						pDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = p[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDot[currlin+j] += migrates;
						pDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = p[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = p[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDot[currlin+j] += migrates;
						pDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
			}		// lineages       		
//		}

		pDot[length-1] /= 2;

	}

void Euler2ndOrderCPU8::computeSecondDerivateWithMultiplicator(double * p, double * pDot, double * pDotDot, int length) {
		//double * sumDotStates = new double[states];;
		memset(sumDotStates, 0, states * sizeof(double));

		int j;
		for (int i = 0; i<lineages; i++) {
			for (j = 0; j<states; j++) {
				sumDotStates[j] += multiplicator[i]*pDot[states*i+j];
			}
		}

		// Calculate the change in the lineage state probabilities for every lineage in every state
		for (int i = 0; i<lineages; i++) {
			double pCoalRate = 0.0;
			int currlin = states*i;
			for (int j = 0; j<states; j++)
			pCoalRate += coalescent_rates[j] * (pDot[currlin+j]* (sumStates[j] - p[currlin+j]) + p[currlin+j]* (sumDotStates[j] - pDot[currlin+j]));

			for (j = 0; j<states; j++) {
				pDotDot[currlin+j] *= pDot[currlin+j];
			}

			pDotDot[length-1] -= multiplicator[i]*pCoalRate;

			for (j = 0; j<states; j++) {
				pDotDot[currlin+j] += p[currlin+j]*(pCoalRate - coalescent_rates[j] * (sumDotStates[j] - pDot[currlin+j]));
			}		// j    		
		}		// lineages 

		double migrates;

				// Calculate the probability of a lineage changing states
//		if (hasIndicators) {
//			for (int j = 0; j < indicators.length/2; j++) {
//				int source = indicators[j * n2 + 0];
//				int sink = indicators[j * n2 + 1];
//				double mrate = migration_rates[source * n + sink];
//				for (int i = 0; i<lineages; i++) {
//					migrates = pDot[states*i+source]*mrate;
//					pDotDot[states*i+sink] += migrates;
//					pDotDot[states*i+source] -= migrates;
//				}
//			}
//		} else {
			for (int i = 0; i<lineages; i++) {
				int currlin = states*i;
				// Calculate the probability of a lineage changing states
j = 0;
{
					double pj = pDot[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[currlin+j] += migrates;
						pDotDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = pDot[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[currlin+j] += migrates;
						pDotDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = pDot[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[currlin+j] += migrates;
						pDotDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = pDot[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[currlin+j] += migrates;
						pDotDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = pDot[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[currlin+j] += migrates;
						pDotDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = pDot[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[currlin+j] += migrates;
						pDotDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = pDot[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[currlin+j] += migrates;
						pDotDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
j++;
{
					double pj = pDot[currlin+j];
					for (int k = j+1; k < states; k++) {

						// the probability of lineage i being in state j is p[i*nr_states +j]
						migrates = pDot[currlin+k]*migration_rates[k * n + j] -
						pj*migration_rates[j * n + k];
						pDotDot[currlin+j] += migrates;
						pDotDot[currlin+k] -= migrates;
					}		// j XXX

				}		// j
			}		// lineages    

//		}
		pDotDot[length-1] /= 2;
	}



