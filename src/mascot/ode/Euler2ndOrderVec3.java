package mascot.ode;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.apache.commons.math3.util.FastMath;

/**
 * SIMD version of Euler2ndOrderFused3 (3 states, canonical indicator order) using the Vector API:
 * several lineages are processed per instruction. Lineage state probabilities are converted to
 * one array per state on entry and back on exit. Sums over lineages are accumulated per lane, so
 * results agree with Euler2ndOrderFused3 up to floating point rounding, but are not bit-identical.
 * Requires the JVM option --add-modules jdk.incubator.vector.
 */
public class Euler2ndOrderVec3 extends Euler2ndOrderFused3 {

	static final VectorSpecies<Double> SP = DoubleVector.SPECIES_PREFERRED;
	static final int W = SP.length();

	// one array per state, padded to a multiple of W with zero probability lineages
	double[] P0 = new double[0], P1, P2, D0, D1, D2, DD0, DD1, DD2;

	private void ensureCapacity(int n) {
		if (P0.length >= n) return;
		P0 = new double[n]; P1 = new double[n]; P2 = new double[n];
		D0 = new double[n]; D1 = new double[n]; D2 = new double[n];
		DD0 = new double[n]; DD1 = new double[n]; DD2 = new double[n];
	}

	@Override
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
		final int L = lineages;
		if (hasMultiplicator || states != 3 || !hasIndicators || L < 2 * W || !canonicalRates()) {
			super.calculateValues(duration, p, pDot, pDotDot, pDotDotDot, length);
			return;
		}
		final int Lpad = ((L + W - 1) / W) * W;
		ensureCapacity(Lpad);
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2;
		final double[] D0 = this.D0, D1 = this.D1, D2 = this.D2;
		final double[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2;
		for (int i = 0, k = 0; i < L; i++, k += 3) {
			P0[i] = p[k]; P1[i] = p[k + 1]; P2[i] = p[k + 2];
		}
		for (int i = L; i < Lpad; i++) {
			P0[i] = 0; P1[i] = 0; P2[i] = 0;
		}
		final int last = length - 1;
		double logPSlot = p[last];

		final double cr0 = coalescent_rates[0], cr1 = coalescent_rates[1], cr2 = coalescent_rates[2];
		final double r01 = r[0], r02 = r[1], r10 = r[2], r12 = r[3], r20 = r[4], r21 = r[5];
		final DoubleVector vcr0 = DoubleVector.broadcast(SP, cr0), vcr1 = DoubleVector.broadcast(SP, cr1), vcr2 = DoubleVector.broadcast(SP, cr2);
		final DoubleVector vr01 = DoubleVector.broadcast(SP, r01), vr02 = DoubleVector.broadcast(SP, r02);
		final DoubleVector vr10 = DoubleVector.broadcast(SP, r10), vr12 = DoubleVector.broadcast(SP, r12);
		final DoubleVector vr20 = DoubleVector.broadcast(SP, r20), vr21 = DoubleVector.broadcast(SP, r21);
		// total rate out of each state
		final DoubleVector vout0 = DoubleVector.broadcast(SP, r01 + r02);
		final DoubleVector vout1 = DoubleVector.broadcast(SP, r10 + r12);
		final DoubleVector vout2 = DoubleVector.broadcast(SP, r20 + r21);
		final DoubleVector zero = DoubleVector.zero(SP);
		final int lastBlock = Lpad - W;
		final VectorMask<Double> lastValid = SP.indexInRange(lastBlock, L);

		double sS0 = 0, sS1 = 0, sS2 = 0;
		for (int i = 0; i < L; i++) {
			sS0 += P0[i]; sS1 += P1[i]; sS2 += P2[i];
		}

		while (duration > 0) {
			iterations++;
			final DoubleVector vsS0 = DoubleVector.broadcast(SP, sS0), vsS1 = DoubleVector.broadcast(SP, sS1), vsS2 = DoubleVector.broadcast(SP, sS2);

			// pass 1: first derivative, sums of first derivatives
			DoubleVector accCoal = zero, accD0 = zero, accD1 = zero, accD2 = zero;
			for (int i = 0; i < Lpad; i += W) {
				DoubleVector p0 = DoubleVector.fromArray(SP, P0, i);
				DoubleVector p1 = DoubleVector.fromArray(SP, P1, i);
				DoubleVector p2 = DoubleVector.fromArray(SP, P2, i);
				DoubleVector t0 = vsS0.sub(p0).mul(vcr0), t1 = vsS1.sub(p1).mul(vcr1), t2 = vsS2.sub(p2).mul(vcr2);
				DoubleVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2));
				accCoal = accCoal.add(sc);
				DoubleVector d0 = p0.mul(sc.sub(t0)).add(p1.mul(vr10)).add(p2.mul(vr20)).sub(p0.mul(vout0));
				DoubleVector d1 = p1.mul(sc.sub(t1)).add(p0.mul(vr01)).add(p2.mul(vr21)).sub(p1.mul(vout1));
				DoubleVector d2 = p2.mul(sc.sub(t2)).add(p0.mul(vr02)).add(p1.mul(vr12)).sub(p2.mul(vout2));
				d0.intoArray(D0, i); d1.intoArray(D1, i); d2.intoArray(D2, i);
				accD0 = accD0.add(d0); accD1 = accD1.add(d1); accD2 = accD2.add(d2);
			}
			final double dLast = -accCoal.reduceLanes(VectorOperators.ADD) / 2;
			final double sD0 = accD0.reduceLanes(VectorOperators.ADD), sD1 = accD1.reduceLanes(VectorOperators.ADD), sD2 = accD2.reduceLanes(VectorOperators.ADD);
			final DoubleVector vsD0 = DoubleVector.broadcast(SP, sD0), vsD1 = DoubleVector.broadcast(SP, sD1), vsD2 = DoubleVector.broadcast(SP, sD2);

			// pass 2: second and approximate third derivative, max of the third derivative
			DoubleVector accPc = zero, accMax = zero;
			for (int i = 0; i < Lpad; i += W) {
				DoubleVector p0 = DoubleVector.fromArray(SP, P0, i);
				DoubleVector p1 = DoubleVector.fromArray(SP, P1, i);
				DoubleVector p2 = DoubleVector.fromArray(SP, P2, i);
				DoubleVector d0 = DoubleVector.fromArray(SP, D0, i);
				DoubleVector d1 = DoubleVector.fromArray(SP, D1, i);
				DoubleVector d2 = DoubleVector.fromArray(SP, D2, i);
				DoubleVector t0 = vsS0.sub(p0).mul(vcr0), t1 = vsS1.sub(p1).mul(vcr1), t2 = vsS2.sub(p2).mul(vcr2);
				DoubleVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2));
				DoubleVector c0 = sc.sub(t0), c1 = sc.sub(t1), c2 = sc.sub(t2);

				DoubleVector pc = vcr0.mul(d0.mul(vsS0.sub(p0.mul(2))).add(p0.mul(vsD0)))
						.add(vcr1.mul(d1.mul(vsS1.sub(p1.mul(2))).add(p1.mul(vsD1))))
						.add(vcr2.mul(d2.mul(vsS2.sub(p2.mul(2))).add(p2.mul(vsD2))));
				accPc = accPc.add(pc);
				DoubleVector dd0 = c0.mul(d0).add(p0.mul(pc.sub(vcr0.mul(vsD0.sub(d0)))))
						.add(d1.mul(vr10)).add(d2.mul(vr20)).sub(d0.mul(vout0));
				DoubleVector dd1 = c1.mul(d1).add(p1.mul(pc.sub(vcr1.mul(vsD1.sub(d1)))))
						.add(d0.mul(vr01)).add(d2.mul(vr21)).sub(d1.mul(vout1));
				DoubleVector dd2 = c2.mul(d2).add(p2.mul(pc.sub(vcr2.mul(vsD2.sub(d2)))))
						.add(d0.mul(vr02)).add(d1.mul(vr12)).sub(d2.mul(vout2));
				dd0.intoArray(DD0, i); dd1.intoArray(DD1, i); dd2.intoArray(DD2, i);

				DoubleVector e0 = c0.mul(dd0).add(dd1.mul(vr10)).add(dd2.mul(vr20)).sub(dd0.mul(vout0));
				DoubleVector e1 = c1.mul(dd1).add(dd0.mul(vr01)).add(dd2.mul(vr21)).sub(dd1.mul(vout1));
				DoubleVector e2 = c2.mul(dd2).add(dd0.mul(vr02)).add(dd1.mul(vr12)).sub(dd2.mul(vout2));
				accMax = accMax.max(e0.abs()).max(e1.abs()).max(e2.abs());
			}
			final double ddLast = -accPc.reduceLanes(VectorOperators.ADD) / 2;
			final double maxDDD = accMax.reduceLanes(VectorOperators.MAX);

			// pass 3: time step; check all elements at once, and only search a smaller step if needed
			double timeStep = FastMath.min(FastMath.cbrt(epsilon * 6 / maxDDD), FastMath.min(duration, max_step));
			double timeStepSquare = timeStep * timeStep * 0.5;
			boolean needsSearch = false;
			{
				final DoubleVector vts = DoubleVector.broadcast(SP, timeStep), vtss = DoubleVector.broadcast(SP, timeStepSquare);
				for (int i = 0; i < Lpad && !needsSearch; i += W) {
					needsSearch = outOfRange(P0, D0, DD0, i, vts, vtss) || outOfRange(P1, D1, DD1, i, vts, vtss) || outOfRange(P2, D2, DD2, i, vts, vtss);
				}
			}
			if (needsSearch) {
				// same sequential search as Euler2ndOrder3.updateP, in the original element order
				search:
				for (int i = 0; i < L; i++) {
					for (int j = 0; j < 3; j++) {
						final double[] PJ = j == 0 ? P0 : (j == 1 ? P1 : P2);
						final double[] DJ = j == 0 ? D0 : (j == 1 ? D1 : D2);
						final double[] DDJ = j == 0 ? DD0 : (j == 1 ? DD1 : DD2);
						double new_val = PJ[i] + DJ[i] * timeStep + DDJ[i] * timeStepSquare;
						double diff = FastMath.abs(new_val - PJ[i]);
						int its = 0;
						while (new_val > 1 || new_val < 0 || diff > 0.2) {
							timeStep *= 0.9;
							timeStepSquare = timeStep * timeStep * 0.5;
							new_val = PJ[i] + DJ[i] * timeStep + DDJ[i] * timeStepSquare;
							diff = FastMath.abs(new_val - PJ[i]);
							its++;
							if (its > 10000) {
								P2[L - 1] = Double.NEGATIVE_INFINITY;
								break search;
							}
						}
					}
				}
			}

			if (P2[L - 1] == Double.NEGATIVE_INFINITY) {
				duration = 0.0;
			} else {
				// pass 4: update, normalise, sum states for the next step
				logPSlot += dLast * timeStep + ddLast * timeStepSquare;
				final DoubleVector vts = DoubleVector.broadcast(SP, timeStep), vtss = DoubleVector.broadcast(SP, timeStepSquare);
				DoubleVector accS0 = zero, accS1 = zero, accS2 = zero;
				for (int i = 0; i < Lpad; i += W) {
					DoubleVector x0 = DoubleVector.fromArray(SP, P0, i).add(DoubleVector.fromArray(SP, D0, i).mul(vts).add(DoubleVector.fromArray(SP, DD0, i).mul(vtss)));
					DoubleVector x1 = DoubleVector.fromArray(SP, P1, i).add(DoubleVector.fromArray(SP, D1, i).mul(vts).add(DoubleVector.fromArray(SP, DD1, i).mul(vtss)));
					DoubleVector x2 = DoubleVector.fromArray(SP, P2, i).add(DoubleVector.fromArray(SP, D2, i).mul(vts).add(DoubleVector.fromArray(SP, DD2, i).mul(vtss)));
					VectorMask<Double> negative = x0.compare(VectorOperators.LT, 0.0).or(x1.compare(VectorOperators.LT, 0.0)).or(x2.compare(VectorOperators.LT, 0.0));
					DoubleVector sum = x0.add(x1).add(x2);
					DoubleVector y0 = x0.div(sum), y1 = x1.div(sum), y2 = x2.div(sum);
					if (negative.anyTrue()) {
						// like normalise(): flag the failure and leave the lineage unnormalised
						logPSlot = Double.NEGATIVE_INFINITY;
						y0 = y0.blend(x0, negative); y1 = y1.blend(x1, negative); y2 = y2.blend(x2, negative);
					}
					if (i == lastBlock) {
						// padding lineages stay at zero (0/0 would give NaN)
						VectorMask<Double> pad = lastValid.not();
						y0 = y0.blend(0.0, pad); y1 = y1.blend(0.0, pad); y2 = y2.blend(0.0, pad);
					}
					y0.intoArray(P0, i); y1.intoArray(P1, i); y2.intoArray(P2, i);
					accS0 = accS0.add(y0); accS1 = accS1.add(y1); accS2 = accS2.add(y2);
				}
				sS0 = accS0.reduceLanes(VectorOperators.ADD);
				sS1 = accS1.reduceLanes(VectorOperators.ADD);
				sS2 = accS2.reduceLanes(VectorOperators.ADD);
				duration -= timeStep;
			}

			if (iterations > 10000) {
				System.err.println("too many iterations, return negative infinity");
				logPSlot = Double.NEGATIVE_INFINITY;
				break;
			}
		}

		for (int i = 0, k = 0; i < L; i++, k += 3) {
			p[k] = P0[i]; p[k + 1] = P1[i]; p[k + 2] = P2[i];
		}
		p[last] = logPSlot;
	}

	private static boolean outOfRange(double[] P, double[] D, double[] DD, int i, DoubleVector vts, DoubleVector vtss) {
		DoubleVector p = DoubleVector.fromArray(SP, P, i);
		DoubleVector n = p.add(DoubleVector.fromArray(SP, D, i).mul(vts)).add(DoubleVector.fromArray(SP, DD, i).mul(vtss));
		return n.compare(VectorOperators.GT, 1.0)
				.or(n.compare(VectorOperators.LT, 0.0))
				.or(n.sub(p).abs().compare(VectorOperators.GT, 0.2))
				.anyTrue();
	}
}
