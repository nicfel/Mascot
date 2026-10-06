package mascot.ode;

import org.apache.commons.math3.util.FastMath;

/**
 * Loop-fused version of Euler2ndOrder: performs the same floating point operations
 * in the same order per element (so results are bit-identical to Euler2ndOrder/Euler2ndOrderN),
 * but touches memory 4 times per integration step instead of ~20 times, and works lineage by lineage.
 * The multiplicator path is delegated to the original implementation.
 */
public class Euler2ndOrderFused extends Euler2ndOrder {

	// migration terms of the current rates interval, in the order the original code applies them
	int[] pairSrc = new int[0], pairSnk = new int[0];
	double[] pairRate = new double[0];
	int nPairs;
	double[] sumS, sumD, tcr, pk;

	@Override
	public void setup(int maxSize, int states, double epsilon, double max_step) {
		super.setup(maxSize, states, epsilon, max_step);
		sumS = new double[states];
		sumD = new double[states];
		tcr = new double[states];
	}

	@Override
	public void initAndcalculateValues(int ratesInterval, int lineages, double duration, double[] p, int length) {
		if (ratesInterval >= nextRateShift.length) {
			ratesInterval = nextRateShift.length - 1;
		}
		migration_rates = migrationRates[ratesInterval];
		coalescent_rates = coalescentRates[ratesInterval];
		indicators = indicators_[ratesInterval];
		iterations = 0;
		n = (int) (Math.sqrt(migration_rates.length) + 0.5);
		this.lineages = lineages;
		this.dimension = lineages * states;
		hasIndicators = (indicators != null);
		hasMultiplicator = false;
		calculateValues(duration, p, linProbs_tmpdt, linProbs_tmpddt, linProbs_tmpdddt, length);
	}

	private void setUpPairs() {
		int max = hasIndicators ? indicators.length / 2 : states * (states - 1) / 2;
		if (pairSrc.length < max) {
			pairSrc = new int[max]; pairSnk = new int[max]; pairRate = new double[max];
		}
		nPairs = 0;
		if (hasIndicators) {
			for (int j = 0; j < indicators.length / 2; j++) {
				int source = indicators[j * n2 + 0];
				int sink = indicators[j * n2 + 1];
				pairSrc[nPairs] = source; pairSnk[nPairs] = sink;
				pairRate[nPairs] = migration_rates[source * n + sink];
				nPairs++;
			}
		}
	}

	@Override
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
		if (hasMultiplicator) {
			super.calculateValues(duration, p, pDot, pDotDot, pDotDotDot, length);
			return;
		}
		setUpPairs();
		final int S = states, L = lineages, last = length - 1;
		final double[] cr = coalescent_rates, M = migration_rates;
		final double[] sumS = this.sumS, sumD = this.sumD, tcr = this.tcr;

		for (int j = 0; j < S; j++) sumS[j] = 0.0;
		for (int k = 0; k < L * S; k += S) {
			for (int j = 0; j < S; j++) sumS[j] += p[k + j];
		}

		while (duration > 0) {
			iterations++;

			// pass 1: first derivative per lineage, sum of first derivatives
			double dLast = 0.0;
			for (int j = 0; j < S; j++) sumD[j] = 0.0;
			for (int k = 0; k < L * S; k += S) {
				double sumCoal = 0;
				for (int j = 0; j < S; j++) {
					tcr[j] = cr[j] * (sumS[j] - p[k + j]);
					sumCoal += p[k + j] * tcr[j];
				}
				dLast -= sumCoal;
				for (int j = 0; j < S; j++) {
					double coal = sumCoal - tcr[j];
					pDotDotDot[k + j] = coal;
					pDot[k + j] = 0.0 + p[k + j] * coal;
				}
				migrate(p, pDot, k, M);
				for (int j = 0; j < S; j++) sumD[j] += pDot[k + j];
			}
			pDot[last] = dLast / 2;

			// pass 2: second and approximate third derivative per lineage, max of the third derivative
			double ddLast = 0.0;
			double maxDDD = 0.0;
			for (int k = 0; k < L * S; k += S) {
				double pCoalRate = 0.0;
				for (int j = 0; j < S; j++) {
					pCoalRate += cr[j] * (pDot[k + j] * (sumS[j] - 2 * p[k + j]) + p[k + j] * (sumD[j]));
				}
				for (int j = 0; j < S; j++) {
					pDotDot[k + j] = pDotDotDot[k + j] * pDot[k + j] + p[k + j] * (pCoalRate - cr[j] * (sumD[j] - pDot[k + j]));
				}
				ddLast -= pCoalRate;
				migrate(pDot, pDotDot, k, M);
				for (int j = 0; j < S; j++) pDotDotDot[k + j] *= pDotDot[k + j];
				migrateThird(pDotDot, pDotDotDot, k, M);
				for (int j = 0; j < S; j++) maxDDD = FastMath.max(maxDDD, FastMath.abs(pDotDotDot[k + j]));
			}
			pDotDot[last] = ddLast / 2;

			// pass 3: find the time step (identical to updateP)
			final int len = last;
			double timeStep = FastMath.min(FastMath.cbrt(epsilon * 6 / maxDDD), FastMath.min(duration, max_step));
			double timeStepSquare = timeStep * timeStep * 0.5;
			for (int i = 0; i < len; i++) {
				double new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
				double diff = FastMath.abs(new_val - p[i]);
				int its = 0;
				while (new_val > 1 || new_val < 0 || diff > 0.2) {
					timeStep *= 0.9;
					timeStepSquare = timeStep * timeStep * 0.5;
					new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
					diff = FastMath.abs(new_val - p[i]);
					its++;
					if (its > 10000) {
						p[len - 1] = Double.NEGATIVE_INFINITY;
						break;
					}
				}
			}
			if (p[len - 1] == Double.NEGATIVE_INFINITY) {
				duration = 0.0;
			} else {
				// pass 4: update, normalise, and sum states for the next step
				p[last] += pDot[last] * timeStep + pDotDot[last] * timeStepSquare;
				for (int j = 0; j < S; j++) sumS[j] = 0.0;
				for (int k = 0; k < L * S; k += S) {
					double linSum = 0;
					boolean negative = false;
					for (int j = 0; j < S; j++) {
						double x = p[k + j] + (pDot[k + j] * timeStep + pDotDot[k + j] * timeStepSquare);
						p[k + j] = x;
					}
					for (int j = 0; j < S; j++) {
						double x = p[k + j];
						linSum += x;
						if (x < 0.0) {
							p[last] = Double.NEGATIVE_INFINITY;
							negative = true;
							break;
						}
					}
					if (!negative) {
						for (int j = 0; j < S; j++) p[k + j] /= linSum;
					}
					for (int j = 0; j < S; j++) sumS[j] += p[k + j];
				}
				duration -= timeStep;
			}

			if (iterations > 10000) {
				System.err.println("too many iterations, return negative infinity");
				p[length - 1] = Double.NEGATIVE_INFINITY;
				break;
			}
		}
	}

	/** apply migration to one lineage starting at k: to[] += terms of from[], in the order of the original code */
	private void migrate(final double[] from, final double[] to, final int k, final double[] M) {
		if (hasIndicators) {
			for (int q = 0; q < nPairs; q++) {
				double migrates = from[k + pairSrc[q]] * pairRate[q];
				to[k + pairSnk[q]] += migrates;
				to[k + pairSrc[q]] -= migrates;
			}
		} else {
			for (int j = 0; j < states; j++) {
				double pj = from[k + j];
				for (int l = j + 1; l < states; l++) {
					double migrates = from[k + l] * M[l * n + j] - pj * M[j * n + l];
					to[k + j] += migrates;
					to[k + l] -= migrates;
				}
			}
		}
	}

	private void migrateThird(final double[] pDotDot, final double[] pDotDotDot, final int k, final double[] M) {
		if (hasIndicators) {
			for (int q = 0; q < nPairs; q++) {
				double migrates = pDotDot[k + pairSrc[q]] * pairRate[q];
				pDotDotDot[k + pairSnk[q]] += migrates;
				pDotDotDot[k + pairSrc[q]] -= migrates;
			}
		} else {
			for (int j = 0; j < states; j++) {
				for (int l = 0; l < states; l++) {
					double migrates = pDotDot[k + j] * M[j * n + l];
					pDotDotDot[k + l] += migrates;
					pDotDotDot[k + j] -= migrates;
				}
			}
		}
	}
}
