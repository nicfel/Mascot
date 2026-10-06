package mascot.ode;

import org.apache.commons.math3.util.FastMath;

/**
 * Unrolled 3-state version of Euler2ndOrderFused. Same per-element operation order as Euler2ndOrder3
 * (bit-identical results). Indicator pairs must come in the canonical order
 * (0,1),(0,2),(1,0),(1,2),(2,0),(2,1) (as produced by the dynamics); missing pairs get rate 0,
 * otherwise this falls back to the generic fused code.
 */
public class Euler2ndOrderFused3 extends Euler2ndOrderFused {

	static final int[] CANON_SRC = {0, 0, 1, 1, 2, 2};
	static final int[] CANON_SNK = {1, 2, 0, 2, 0, 1};
	final double[] r = new double[6];

	/** fills r with canonical pair rates; returns false if the indicator order is not canonical */
	protected boolean canonicalRates() {
		java.util.Arrays.fill(r, 0.0);
		int slot = 0;
		for (int q = 0; q < indicators.length / 2; q++) {
			int src = indicators[q * n2], snk = indicators[q * n2 + 1];
			if (src == snk) {
				if (migration_rates[src * n + snk] != 0.0) return false;
				continue;
			}
			while (slot < 6 && (CANON_SRC[slot] != src || CANON_SNK[slot] != snk)) slot++;
			if (slot == 6) return false;
			r[slot++] = migration_rates[src * n + snk];
		}
		return true;
	}

	@Override
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
		if (hasMultiplicator || states != 3 || !hasIndicators || !canonicalRates()) {
			super.calculateValues(duration, p, pDot, pDotDot, pDotDotDot, length);
			return;
		}
		final int L3 = lineages * 3, last = length - 1;
		final double cr0 = coalescent_rates[0], cr1 = coalescent_rates[1], cr2 = coalescent_rates[2];
		final double r01 = r[0], r02 = r[1], r10 = r[2], r12 = r[3], r20 = r[4], r21 = r[5];

		double sS0 = 0.0, sS1 = 0.0, sS2 = 0.0;
		for (int k = 0; k < L3; k += 3) {
			sS0 += p[k]; sS1 += p[k + 1]; sS2 += p[k + 2];
		}

		while (duration > 0) {
			iterations++;

			// pass 1: first derivative, sums of first derivative
			double dLast = 0.0, sD0 = 0.0, sD1 = 0.0, sD2 = 0.0;
			for (int k = 0; k < L3; k += 3) {
				final double p0 = p[k], p1 = p[k + 1], p2 = p[k + 2];
				final double t0 = cr0 * (sS0 - p0), t1 = cr1 * (sS1 - p1), t2 = cr2 * (sS2 - p2);
				double sumCoal = 0;
				sumCoal += p0 * t0; sumCoal += p1 * t1; sumCoal += p2 * t2;
				dLast -= sumCoal;
				double d0 = 0.0 + p0 * (sumCoal - t0);
				double d1 = 0.0 + p1 * (sumCoal - t1);
				double d2 = 0.0 + p2 * (sumCoal - t2);
				double m;
				m = p0 * r01; d1 += m; d0 -= m;
				m = p0 * r02; d2 += m; d0 -= m;
				m = p1 * r10; d0 += m; d1 -= m;
				m = p1 * r12; d2 += m; d1 -= m;
				m = p2 * r20; d0 += m; d2 -= m;
				m = p2 * r21; d1 += m; d2 -= m;
				pDot[k] = d0; pDot[k + 1] = d1; pDot[k + 2] = d2;
				sD0 += d0; sD1 += d1; sD2 += d2;
			}
			pDot[last] = dLast / 2;

			// pass 2: second derivative, approximate third derivative and its max
			double ddLast = 0.0;
			double maxDDD = 0.0;
			for (int k = 0; k < L3; k += 3) {
				final double p0 = p[k], p1 = p[k + 1], p2 = p[k + 2];
				final double d0 = pDot[k], d1 = pDot[k + 1], d2 = pDot[k + 2];
				// recompute the coalescent terms of pass 1 (same inputs, same result)
				final double t0 = cr0 * (sS0 - p0), t1 = cr1 * (sS1 - p1), t2 = cr2 * (sS2 - p2);
				double sumCoal = 0;
				sumCoal += p0 * t0; sumCoal += p1 * t1; sumCoal += p2 * t2;
				final double c0 = sumCoal - t0, c1 = sumCoal - t1, c2 = sumCoal - t2;

				double pc = 0.0;
				pc += cr0 * (d0 * (sS0 - 2 * p0) + p0 * (sD0));
				pc += cr1 * (d1 * (sS1 - 2 * p1) + p1 * (sD1));
				pc += cr2 * (d2 * (sS2 - 2 * p2) + p2 * (sD2));
				double dd0 = c0 * d0 + p0 * (pc - cr0 * (sD0 - d0));
				double dd1 = c1 * d1 + p1 * (pc - cr1 * (sD1 - d1));
				double dd2 = c2 * d2 + p2 * (pc - cr2 * (sD2 - d2));
				ddLast -= pc;
				double m;
				m = d0 * r01; dd1 += m; dd0 -= m;
				m = d0 * r02; dd2 += m; dd0 -= m;
				m = d1 * r10; dd0 += m; dd1 -= m;
				m = d1 * r12; dd2 += m; dd1 -= m;
				m = d2 * r20; dd0 += m; dd2 -= m;
				m = d2 * r21; dd1 += m; dd2 -= m;
				pDotDot[k] = dd0; pDotDot[k + 1] = dd1; pDotDot[k + 2] = dd2;

				double e0 = c0 * dd0, e1 = c1 * dd1, e2 = c2 * dd2;
				m = dd0 * r01; e1 += m; e0 -= m;
				m = dd0 * r02; e2 += m; e0 -= m;
				m = dd1 * r10; e0 += m; e1 -= m;
				m = dd1 * r12; e2 += m; e1 -= m;
				m = dd2 * r20; e0 += m; e2 -= m;
				m = dd2 * r21; e1 += m; e2 -= m;
				maxDDD = Math.max(maxDDD, Math.abs(e0));
				maxDDD = Math.max(maxDDD, Math.abs(e1));
				maxDDD = Math.max(maxDDD, Math.abs(e2));
			}
			pDotDot[last] = ddLast / 2;

			// pass 3: find the time step (identical to updateP)
			double timeStep = FastMath.min(FastMath.cbrt(epsilon * 6 / maxDDD), FastMath.min(duration, max_step));
			double timeStepSquare = timeStep * timeStep * 0.5;
			for (int i = 0; i < L3; i++) {
				double new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
				if (new_val > 1 || new_val < 0 || Math.abs(new_val - p[i]) > 0.2) {
					double diff = Math.abs(new_val - p[i]);
					int its = 0;
					while (new_val > 1 || new_val < 0 || diff > 0.2) {
						timeStep *= 0.9;
						timeStepSquare = timeStep * timeStep * 0.5;
						new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
						diff = Math.abs(new_val - p[i]);
						its++;
						if (its > 10000) {
							p[L3 - 1] = Double.NEGATIVE_INFINITY;
							break;
						}
					}
				}
			}
			if (p[L3 - 1] == Double.NEGATIVE_INFINITY) {
				duration = 0.0;
			} else {
				// pass 4: update, normalise, sum states for the next step
				p[last] += pDot[last] * timeStep + pDotDot[last] * timeStepSquare;
				sS0 = 0.0; sS1 = 0.0; sS2 = 0.0;
				for (int k = 0; k < L3; k += 3) {
					double x0 = p[k] + (pDot[k] * timeStep + pDotDot[k] * timeStepSquare);
					double x1 = p[k + 1] + (pDot[k + 1] * timeStep + pDotDot[k + 1] * timeStepSquare);
					double x2 = p[k + 2] + (pDot[k + 2] * timeStep + pDotDot[k + 2] * timeStepSquare);
					if (x0 < 0.0 || x1 < 0.0 || x2 < 0.0) {
						// original normalise() flags the failure and leaves the lineage unnormalised
						p[last] = Double.NEGATIVE_INFINITY;
					} else {
						double linSum = 0;
						linSum += x0; linSum += x1; linSum += x2;
						x0 /= linSum; x1 /= linSum; x2 /= linSum;
					}
					p[k] = x0; p[k + 1] = x1; p[k + 2] = x2;
					sS0 += x0; sS1 += x1; sS2 += x2;
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
}
