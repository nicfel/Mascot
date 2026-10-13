package mascot.ode;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.apache.commons.math3.util.FastMath;

/**
 * SIMD (Vector API) solver for 5 states in double precision, written out like Euler2ndOrderVec3:
 * several lineages are processed per instruction, per-state values stay in vector registers and
 * all migration terms are unrolled. Lineage state probabilities are kept in one array per state;
 * sums over lineages are accumulated per lane, so results agree with Euler2ndOrder5 up to rounding.
 * Indicator pairs must come in the canonical order (a, b != a) as produced by the dynamics,
 * otherwise (and for very few lineages) Euler2ndOrderFused is used.
 * Requires the JVM option --add-modules jdk.incubator.vector.
 *
 * Generated code, do not edit by hand.
 */
public class Euler2ndOrderVec5 extends Euler2ndOrderFused {

	static final int N = 5;
	static final int NP = N * (N - 1);
	static final VectorSpecies<Double> SP = DoubleVector.SPECIES_PREFERRED;
	static final int W = SP.length();
	// migration routes (source, sink) in canonical order
	static final int[] SRC = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
	static final int[] SNK = {1, 2, 3, 4, 0, 2, 3, 4, 0, 1, 3, 4, 0, 1, 2, 4, 0, 1, 2, 3};

	final double[] rt = new double[NP], out = new double[N], cr = new double[N];
	final double[] sS = new double[N], sD = new double[N];
	double dLast, ddLast, maxDDD;
	double[] P0 = new double[0], P1 = new double[0], P2 = new double[0], P3 = new double[0], P4 = new double[0], D0, D1, D2, D3, D4, DD0, DD1, DD2, DD3, DD4;

	private void ensureCapacity(int size) {
		if (P0.length >= size) return;
		P0 = new double[size]; D0 = new double[size]; DD0 = new double[size];
		P1 = new double[size]; D1 = new double[size]; DD1 = new double[size];
		P2 = new double[size]; D2 = new double[size]; DD2 = new double[size];
		P3 = new double[size]; D3 = new double[size]; DD3 = new double[size];
		P4 = new double[size]; D4 = new double[size]; DD4 = new double[size];
	}

	/** fills rt with the migration rates in canonical order; returns false if the indicator order is not canonical */
	private boolean canonicalRates() {
		java.util.Arrays.fill(rt, (double) 0);
		int slot = 0;
		for (int k = 0; k < indicators.length / 2; k++) {
			int src = indicators[k * n2], snk = indicators[k * n2 + 1];
			if (src == snk) {
				if (migration_rates[src * n + snk] != 0.0) return false;
				continue;
			}
			while (slot < NP && (SRC[slot] != src || SNK[slot] != snk)) slot++;
			if (slot == NP) return false;
			rt[slot++] = (double) migration_rates[src * n + snk];
		}
		java.util.Arrays.fill(out, (double) 0);
		for (int k = 0; k < NP; k++) out[SRC[k]] += rt[k];
		return true;
	}

	/** first derivative and its sums over lineages */
	private void pass1(final int Lpad) {
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4;
		final double[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4;
		final double[] rt = this.rt, out = this.out;
		final double s0 = (double) sS[0], s1 = (double) sS[1], s2 = (double) sS[2], s3 = (double) sS[3], s4 = (double) sS[4];
		final double c0 = cr[0], c1 = cr[1], c2 = cr[2], c3 = cr[3], c4 = cr[4];
		DoubleVector accCoal = DoubleVector.zero(SP);
		DoubleVector aD0 = DoubleVector.zero(SP), aD1 = DoubleVector.zero(SP), aD2 = DoubleVector.zero(SP), aD3 = DoubleVector.zero(SP), aD4 = DoubleVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final DoubleVector p0 = DoubleVector.fromArray(SP, P0, i);
			final DoubleVector p1 = DoubleVector.fromArray(SP, P1, i);
			final DoubleVector p2 = DoubleVector.fromArray(SP, P2, i);
			final DoubleVector p3 = DoubleVector.fromArray(SP, P3, i);
			final DoubleVector p4 = DoubleVector.fromArray(SP, P4, i);
			final DoubleVector t0 = p0.neg().add(s0).mul(c0);
			final DoubleVector t1 = p1.neg().add(s1).mul(c1);
			final DoubleVector t2 = p2.neg().add(s2).mul(c2);
			final DoubleVector t3 = p3.neg().add(s3).mul(c3);
			final DoubleVector t4 = p4.neg().add(s4).mul(c4);
			final DoubleVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2)).add(p3.mul(t3)).add(p4.mul(t4));
			accCoal = accCoal.add(sc);
			final DoubleVector d0 = p0.mul(sc.sub(t0)).add(p1.mul(rt[4])).add(p2.mul(rt[8])).add(p3.mul(rt[12])).add(p4.mul(rt[16])).sub(p0.mul(out[0]));
			d0.intoArray(D0, i); aD0 = aD0.add(d0);
			final DoubleVector d1 = p1.mul(sc.sub(t1)).add(p0.mul(rt[0])).add(p2.mul(rt[9])).add(p3.mul(rt[13])).add(p4.mul(rt[17])).sub(p1.mul(out[1]));
			d1.intoArray(D1, i); aD1 = aD1.add(d1);
			final DoubleVector d2 = p2.mul(sc.sub(t2)).add(p0.mul(rt[1])).add(p1.mul(rt[5])).add(p3.mul(rt[14])).add(p4.mul(rt[18])).sub(p2.mul(out[2]));
			d2.intoArray(D2, i); aD2 = aD2.add(d2);
			final DoubleVector d3 = p3.mul(sc.sub(t3)).add(p0.mul(rt[2])).add(p1.mul(rt[6])).add(p2.mul(rt[10])).add(p4.mul(rt[19])).sub(p3.mul(out[3]));
			d3.intoArray(D3, i); aD3 = aD3.add(d3);
			final DoubleVector d4 = p4.mul(sc.sub(t4)).add(p0.mul(rt[3])).add(p1.mul(rt[7])).add(p2.mul(rt[11])).add(p3.mul(rt[15])).sub(p4.mul(out[4]));
			d4.intoArray(D4, i); aD4 = aD4.add(d4);
		}
		dLast = -accCoal.reduceLanes(VectorOperators.ADD) / 2;
		sD[0] = aD0.reduceLanes(VectorOperators.ADD);
		sD[1] = aD1.reduceLanes(VectorOperators.ADD);
		sD[2] = aD2.reduceLanes(VectorOperators.ADD);
		sD[3] = aD3.reduceLanes(VectorOperators.ADD);
		sD[4] = aD4.reduceLanes(VectorOperators.ADD);
	}

	/** second derivative and the max of the approximate third derivative */
	private void pass2(final int Lpad) {
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4;
		final double[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4;
		final double[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4;
		final double[] rt = this.rt, out = this.out;
		final double s0 = (double) sS[0], s1 = (double) sS[1], s2 = (double) sS[2], s3 = (double) sS[3], s4 = (double) sS[4];
		final double c0 = cr[0], c1 = cr[1], c2 = cr[2], c3 = cr[3], c4 = cr[4];
		final double q0 = (double) sD[0], q1 = (double) sD[1], q2 = (double) sD[2], q3 = (double) sD[3], q4 = (double) sD[4];
		DoubleVector accPc = DoubleVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final DoubleVector p0 = DoubleVector.fromArray(SP, P0, i);
			final DoubleVector p1 = DoubleVector.fromArray(SP, P1, i);
			final DoubleVector p2 = DoubleVector.fromArray(SP, P2, i);
			final DoubleVector p3 = DoubleVector.fromArray(SP, P3, i);
			final DoubleVector p4 = DoubleVector.fromArray(SP, P4, i);
			final DoubleVector d0 = DoubleVector.fromArray(SP, D0, i);
			final DoubleVector d1 = DoubleVector.fromArray(SP, D1, i);
			final DoubleVector d2 = DoubleVector.fromArray(SP, D2, i);
			final DoubleVector d3 = DoubleVector.fromArray(SP, D3, i);
			final DoubleVector d4 = DoubleVector.fromArray(SP, D4, i);
			final DoubleVector t0 = p0.neg().add(s0).mul(c0);
			final DoubleVector t1 = p1.neg().add(s1).mul(c1);
			final DoubleVector t2 = p2.neg().add(s2).mul(c2);
			final DoubleVector t3 = p3.neg().add(s3).mul(c3);
			final DoubleVector t4 = p4.neg().add(s4).mul(c4);
			final DoubleVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2)).add(p3.mul(t3)).add(p4.mul(t4));
			final DoubleVector pc = d0.mul(p0.mul(-2.0).add(s0)).add(p0.mul(q0)).mul(c0).add(d1.mul(p1.mul(-2.0).add(s1)).add(p1.mul(q1)).mul(c1)).add(d2.mul(p2.mul(-2.0).add(s2)).add(p2.mul(q2)).mul(c2)).add(d3.mul(p3.mul(-2.0).add(s3)).add(p3.mul(q3)).mul(c3)).add(d4.mul(p4.mul(-2.0).add(s4)).add(p4.mul(q4)).mul(c4));
			accPc = accPc.add(pc);
			final DoubleVector dd0 = sc.sub(t0).mul(d0).add(p0.mul(pc.sub(d0.neg().add(q0).mul(c0)))).add(d1.mul(rt[4])).add(d2.mul(rt[8])).add(d3.mul(rt[12])).add(d4.mul(rt[16])).sub(d0.mul(out[0]));
			dd0.intoArray(DD0, i);
			final DoubleVector dd1 = sc.sub(t1).mul(d1).add(p1.mul(pc.sub(d1.neg().add(q1).mul(c1)))).add(d0.mul(rt[0])).add(d2.mul(rt[9])).add(d3.mul(rt[13])).add(d4.mul(rt[17])).sub(d1.mul(out[1]));
			dd1.intoArray(DD1, i);
			final DoubleVector dd2 = sc.sub(t2).mul(d2).add(p2.mul(pc.sub(d2.neg().add(q2).mul(c2)))).add(d0.mul(rt[1])).add(d1.mul(rt[5])).add(d3.mul(rt[14])).add(d4.mul(rt[18])).sub(d2.mul(out[2]));
			dd2.intoArray(DD2, i);
			final DoubleVector dd3 = sc.sub(t3).mul(d3).add(p3.mul(pc.sub(d3.neg().add(q3).mul(c3)))).add(d0.mul(rt[2])).add(d1.mul(rt[6])).add(d2.mul(rt[10])).add(d4.mul(rt[19])).sub(d3.mul(out[3]));
			dd3.intoArray(DD3, i);
			final DoubleVector dd4 = sc.sub(t4).mul(d4).add(p4.mul(pc.sub(d4.neg().add(q4).mul(c4)))).add(d0.mul(rt[3])).add(d1.mul(rt[7])).add(d2.mul(rt[11])).add(d3.mul(rt[15])).sub(d4.mul(out[4]));
			dd4.intoArray(DD4, i);
		}
		ddLast = -accPc.reduceLanes(VectorOperators.ADD) / 2;
	}

	/** max of the approximate third derivative */
	private void pass2b(final int Lpad) {
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4;
		final double[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4;
		final double[] rt = this.rt, out = this.out;
		final double s0 = (double) sS[0], s1 = (double) sS[1], s2 = (double) sS[2], s3 = (double) sS[3], s4 = (double) sS[4];
		final double c0 = cr[0], c1 = cr[1], c2 = cr[2], c3 = cr[3], c4 = cr[4];
		DoubleVector accMax = DoubleVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final DoubleVector p0 = DoubleVector.fromArray(SP, P0, i);
			final DoubleVector p1 = DoubleVector.fromArray(SP, P1, i);
			final DoubleVector p2 = DoubleVector.fromArray(SP, P2, i);
			final DoubleVector p3 = DoubleVector.fromArray(SP, P3, i);
			final DoubleVector p4 = DoubleVector.fromArray(SP, P4, i);
			final DoubleVector dd0 = DoubleVector.fromArray(SP, DD0, i);
			final DoubleVector dd1 = DoubleVector.fromArray(SP, DD1, i);
			final DoubleVector dd2 = DoubleVector.fromArray(SP, DD2, i);
			final DoubleVector dd3 = DoubleVector.fromArray(SP, DD3, i);
			final DoubleVector dd4 = DoubleVector.fromArray(SP, DD4, i);
			final DoubleVector t0 = p0.neg().add(s0).mul(c0);
			final DoubleVector t1 = p1.neg().add(s1).mul(c1);
			final DoubleVector t2 = p2.neg().add(s2).mul(c2);
			final DoubleVector t3 = p3.neg().add(s3).mul(c3);
			final DoubleVector t4 = p4.neg().add(s4).mul(c4);
			final DoubleVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2)).add(p3.mul(t3)).add(p4.mul(t4));
			accMax = accMax.max(sc.sub(t0).mul(dd0).add(dd1.mul(rt[4])).add(dd2.mul(rt[8])).add(dd3.mul(rt[12])).add(dd4.mul(rt[16])).sub(dd0.mul(out[0])).abs());
			accMax = accMax.max(sc.sub(t1).mul(dd1).add(dd0.mul(rt[0])).add(dd2.mul(rt[9])).add(dd3.mul(rt[13])).add(dd4.mul(rt[17])).sub(dd1.mul(out[1])).abs());
			accMax = accMax.max(sc.sub(t2).mul(dd2).add(dd0.mul(rt[1])).add(dd1.mul(rt[5])).add(dd3.mul(rt[14])).add(dd4.mul(rt[18])).sub(dd2.mul(out[2])).abs());
			accMax = accMax.max(sc.sub(t3).mul(dd3).add(dd0.mul(rt[2])).add(dd1.mul(rt[6])).add(dd2.mul(rt[10])).add(dd4.mul(rt[19])).sub(dd3.mul(out[3])).abs());
			accMax = accMax.max(sc.sub(t4).mul(dd4).add(dd0.mul(rt[3])).add(dd1.mul(rt[7])).add(dd2.mul(rt[11])).add(dd3.mul(rt[15])).sub(dd4.mul(out[4])).abs());
		}
		maxDDD = accMax.reduceLanes(VectorOperators.MAX);
	}

	/** true if any probability would leave [0,1] or change by more than 0.2 with this step */
	private boolean needsSearch(final int Lpad, final double ts, final double tss) {
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4;
		final double[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4;
		final double[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4;
		for (int i = 0; i < Lpad; i += W) {
			VectorMask<Double> m = SP.maskAll(false);
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P0, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D0, i).mul(ts)).add(DoubleVector.fromArray(SP, DD0, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P1, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D1, i).mul(ts)).add(DoubleVector.fromArray(SP, DD1, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P2, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D2, i).mul(ts)).add(DoubleVector.fromArray(SP, DD2, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P3, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D3, i).mul(ts)).add(DoubleVector.fromArray(SP, DD3, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P4, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D4, i).mul(ts)).add(DoubleVector.fromArray(SP, DD4, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			if (m.anyTrue()) return true;
		}
		return false;
	}

	/** update and normalise; returns true if a probability became negative */
	private boolean pass4(final int Lpad, final int L, final double ts, final double tss) {
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4;
		final double[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4;
		final double[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4;
		final int lastBlock = Lpad - W;
		final VectorMask<Double> pad = SP.indexInRange(lastBlock, L).not();
		boolean failed = false;
		DoubleVector aS0 = DoubleVector.zero(SP), aS1 = DoubleVector.zero(SP), aS2 = DoubleVector.zero(SP), aS3 = DoubleVector.zero(SP), aS4 = DoubleVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final DoubleVector x0 = DoubleVector.fromArray(SP, P0, i).add(DoubleVector.fromArray(SP, D0, i).mul(ts).add(DoubleVector.fromArray(SP, DD0, i).mul(tss)));
			final DoubleVector x1 = DoubleVector.fromArray(SP, P1, i).add(DoubleVector.fromArray(SP, D1, i).mul(ts).add(DoubleVector.fromArray(SP, DD1, i).mul(tss)));
			final DoubleVector x2 = DoubleVector.fromArray(SP, P2, i).add(DoubleVector.fromArray(SP, D2, i).mul(ts).add(DoubleVector.fromArray(SP, DD2, i).mul(tss)));
			final DoubleVector x3 = DoubleVector.fromArray(SP, P3, i).add(DoubleVector.fromArray(SP, D3, i).mul(ts).add(DoubleVector.fromArray(SP, DD3, i).mul(tss)));
			final DoubleVector x4 = DoubleVector.fromArray(SP, P4, i).add(DoubleVector.fromArray(SP, D4, i).mul(ts).add(DoubleVector.fromArray(SP, DD4, i).mul(tss)));
			final VectorMask<Double> negative = x0.compare(VectorOperators.LT, 0.0).or(x1.compare(VectorOperators.LT, 0.0)).or(x2.compare(VectorOperators.LT, 0.0)).or(x3.compare(VectorOperators.LT, 0.0)).or(x4.compare(VectorOperators.LT, 0.0));
			final DoubleVector sum = x0.add(x1).add(x2).add(x3).add(x4);
			final boolean anyNegative = negative.anyTrue();
			if (anyNegative) failed = true;
			DoubleVector y0 = x0.div(sum);
			DoubleVector y1 = x1.div(sum);
			DoubleVector y2 = x2.div(sum);
			DoubleVector y3 = x3.div(sum);
			DoubleVector y4 = x4.div(sum);
			if (anyNegative) {
				// like normalise(): leave the lineage unnormalised
				y0 = y0.blend(x0, negative);
				y1 = y1.blend(x1, negative);
				y2 = y2.blend(x2, negative);
				y3 = y3.blend(x3, negative);
				y4 = y4.blend(x4, negative);
			}
			if (i == lastBlock) {
				// padding lineages stay at zero (0/0 would give NaN)
				y0 = y0.blend(0.0, pad);
				y1 = y1.blend(0.0, pad);
				y2 = y2.blend(0.0, pad);
				y3 = y3.blend(0.0, pad);
				y4 = y4.blend(0.0, pad);
			}
			y0.intoArray(P0, i); aS0 = aS0.add(y0);
			y1.intoArray(P1, i); aS1 = aS1.add(y1);
			y2.intoArray(P2, i); aS2 = aS2.add(y2);
			y3.intoArray(P3, i); aS3 = aS3.add(y3);
			y4.intoArray(P4, i); aS4 = aS4.add(y4);
		}
		sS[0] = aS0.reduceLanes(VectorOperators.ADD);
		sS[1] = aS1.reduceLanes(VectorOperators.ADD);
		sS[2] = aS2.reduceLanes(VectorOperators.ADD);
		sS[3] = aS3.reduceLanes(VectorOperators.ADD);
		sS[4] = aS4.reduceLanes(VectorOperators.ADD);
		return failed;
	}

	@Override
	public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
		final int L = lineages;
		if (hasMultiplicator || states != N || !hasIndicators || n != N || L < 2 * W || !canonicalRates()) {
			super.calculateValues(duration, p, pDot, pDotDot, pDotDotDot, length);
			return;
		}
		final int Lpad = ((L + W - 1) / W) * W;
		ensureCapacity(Lpad);
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4;
		for (int i = 0, k = 0; i < L; i++, k += N) {
			P0[i] = (double) p[k + 0];
			P1[i] = (double) p[k + 1];
			P2[i] = (double) p[k + 2];
			P3[i] = (double) p[k + 3];
			P4[i] = (double) p[k + 4];
		}
		for (int i = L; i < Lpad; i++) {
			P0[i] = 0; P1[i] = 0; P2[i] = 0; P3[i] = 0; P4[i] = 0;
		}
		for (int j = 0; j < N; j++) cr[j] = (double) coalescent_rates[j];
		final int last = length - 1;
		double logPSlot = p[last];
		{ double s = 0; for (int i = 0; i < L; i++) s += P0[i]; sS[0] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P1[i]; sS[1] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P2[i]; sS[2] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P3[i]; sS[3] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P4[i]; sS[4] = s; }

		while (duration > 0) {
			iterations++;
			pass1(Lpad);
			pass2(Lpad);
			pass2b(Lpad);
			double timeStep = FastMath.min(FastMath.cbrt(epsilon * 6 / maxDDD), FastMath.min(duration, max_step));
			double timeStepSquare = timeStep * timeStep * 0.5;
			if (needsSearch(Lpad, (double) timeStep, (double) timeStepSquare)) {
				// same sequential search as Euler2ndOrder.updateP, in the original element order
				final double[][] PA = {P0, P1, P2, P3, P4}, DA = {D0, D1, D2, D3, D4}, DDA = {DD0, DD1, DD2, DD3, DD4};
				search:
				for (int i = 0; i < L; i++) {
					for (int j = 0; j < N; j++) {
						double new_val = PA[j][i] + DA[j][i] * timeStep + DDA[j][i] * timeStepSquare;
						double diff = FastMath.abs(new_val - PA[j][i]);
						int its = 0;
						while (new_val > 1 || new_val < 0 || diff > 0.2) {
							timeStep *= 0.9;
							timeStepSquare = timeStep * timeStep * 0.5;
							new_val = PA[j][i] + DA[j][i] * timeStep + DDA[j][i] * timeStepSquare;
							diff = FastMath.abs(new_val - PA[j][i]);
							its++;
							if (its > 10000) {
								P4[L - 1] = Double.NEGATIVE_INFINITY;
								break search;
							}
						}
					}
				}
			}
			if (P4[L - 1] == Double.NEGATIVE_INFINITY) {
				duration = 0.0;
			} else {
				logPSlot += dLast * timeStep + ddLast * timeStepSquare;
				if (pass4(Lpad, L, (double) timeStep, (double) timeStepSquare)) {
					logPSlot = Double.NEGATIVE_INFINITY;
				}
				duration -= timeStep;
			}
			if (iterations > 10000) {
				System.err.println("too many iterations, return negative infinity");
				logPSlot = Double.NEGATIVE_INFINITY;
				break;
			}
		}

		for (int i = 0, k = 0; i < L; i++, k += N) {
			p[k + 0] = P0[i];
			p[k + 1] = P1[i];
			p[k + 2] = P2[i];
			p[k + 3] = P3[i];
			p[k + 4] = P4[i];
		}
		p[last] = logPSlot;
	}
}
