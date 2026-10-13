package mascot.ode;

import jdk.incubator.vector.FloatVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.apache.commons.math3.util.FastMath;

/**
 * SIMD (Vector API) solver for 13 states in single precision, written out like Euler2ndOrderVec3F:
 * several lineages are processed per instruction, per-state values stay in vector registers and
 * all migration terms are unrolled. Lineage state probabilities are kept in one array per state;
 * sums over lineages are accumulated per lane, so results agree with Euler2ndOrder13 up to rounding.
 * Indicator pairs must come in the canonical order (a, b != a) as produced by the dynamics,
 * otherwise (and for very few lineages) Euler2ndOrderFused is used.
 * Requires the JVM option --add-modules jdk.incubator.vector.
 *
 * Generated code, do not edit by hand.
 */
public class Euler2ndOrderVec13F extends Euler2ndOrderFused {

	static final int N = 13;
	static final int NP = N * (N - 1);
	static final VectorSpecies<Float> SP = FloatVector.SPECIES_PREFERRED;
	static final int W = SP.length();
	// migration routes (source, sink) in canonical order
	static final int[] SRC = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12};
	static final int[] SNK = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

	final float[] rt = new float[NP], out = new float[N], cr = new float[N];
	final double[] sS = new double[N], sD = new double[N];
	double dLast, ddLast, maxDDD;
	float[] P0 = new float[0], P1 = new float[0], P2 = new float[0], P3 = new float[0], P4 = new float[0], P5 = new float[0], P6 = new float[0], P7 = new float[0], P8 = new float[0], P9 = new float[0], P10 = new float[0], P11 = new float[0], P12 = new float[0], D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, DD0, DD1, DD2, DD3, DD4, DD5, DD6, DD7, DD8, DD9, DD10, DD11, DD12;

	private void ensureCapacity(int size) {
		if (P0.length >= size) return;
		P0 = new float[size]; D0 = new float[size]; DD0 = new float[size];
		P1 = new float[size]; D1 = new float[size]; DD1 = new float[size];
		P2 = new float[size]; D2 = new float[size]; DD2 = new float[size];
		P3 = new float[size]; D3 = new float[size]; DD3 = new float[size];
		P4 = new float[size]; D4 = new float[size]; DD4 = new float[size];
		P5 = new float[size]; D5 = new float[size]; DD5 = new float[size];
		P6 = new float[size]; D6 = new float[size]; DD6 = new float[size];
		P7 = new float[size]; D7 = new float[size]; DD7 = new float[size];
		P8 = new float[size]; D8 = new float[size]; DD8 = new float[size];
		P9 = new float[size]; D9 = new float[size]; DD9 = new float[size];
		P10 = new float[size]; D10 = new float[size]; DD10 = new float[size];
		P11 = new float[size]; D11 = new float[size]; DD11 = new float[size];
		P12 = new float[size]; D12 = new float[size]; DD12 = new float[size];
	}

	/** fills rt with the migration rates in canonical order; returns false if the indicator order is not canonical */
	private boolean canonicalRates() {
		java.util.Arrays.fill(rt, (float) 0);
		int slot = 0;
		for (int k = 0; k < indicators.length / 2; k++) {
			int src = indicators[k * n2], snk = indicators[k * n2 + 1];
			if (src == snk) {
				if (migration_rates[src * n + snk] != 0.0) return false;
				continue;
			}
			while (slot < NP && (SRC[slot] != src || SNK[slot] != snk)) slot++;
			if (slot == NP) return false;
			rt[slot++] = (float) migration_rates[src * n + snk];
		}
		java.util.Arrays.fill(out, (float) 0);
		for (int k = 0; k < NP; k++) out[SRC[k]] += rt[k];
		return true;
	}

	/** first derivative and its sums over lineages */
	private void pass1(final int Lpad) {
		final float[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12;
		final float[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4, D5 = this.D5, D6 = this.D6, D7 = this.D7, D8 = this.D8, D9 = this.D9, D10 = this.D10, D11 = this.D11, D12 = this.D12;
		final float[] rt = this.rt, out = this.out;
		final float s0 = (float) sS[0], s1 = (float) sS[1], s2 = (float) sS[2], s3 = (float) sS[3], s4 = (float) sS[4], s5 = (float) sS[5], s6 = (float) sS[6], s7 = (float) sS[7], s8 = (float) sS[8], s9 = (float) sS[9], s10 = (float) sS[10], s11 = (float) sS[11], s12 = (float) sS[12];
		final float c0 = cr[0], c1 = cr[1], c2 = cr[2], c3 = cr[3], c4 = cr[4], c5 = cr[5], c6 = cr[6], c7 = cr[7], c8 = cr[8], c9 = cr[9], c10 = cr[10], c11 = cr[11], c12 = cr[12];
		FloatVector accCoal = FloatVector.zero(SP);
		FloatVector aD0 = FloatVector.zero(SP), aD1 = FloatVector.zero(SP), aD2 = FloatVector.zero(SP), aD3 = FloatVector.zero(SP), aD4 = FloatVector.zero(SP), aD5 = FloatVector.zero(SP), aD6 = FloatVector.zero(SP), aD7 = FloatVector.zero(SP), aD8 = FloatVector.zero(SP), aD9 = FloatVector.zero(SP), aD10 = FloatVector.zero(SP), aD11 = FloatVector.zero(SP), aD12 = FloatVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final FloatVector p0 = FloatVector.fromArray(SP, P0, i);
			final FloatVector p1 = FloatVector.fromArray(SP, P1, i);
			final FloatVector p2 = FloatVector.fromArray(SP, P2, i);
			final FloatVector p3 = FloatVector.fromArray(SP, P3, i);
			final FloatVector p4 = FloatVector.fromArray(SP, P4, i);
			final FloatVector p5 = FloatVector.fromArray(SP, P5, i);
			final FloatVector p6 = FloatVector.fromArray(SP, P6, i);
			final FloatVector p7 = FloatVector.fromArray(SP, P7, i);
			final FloatVector p8 = FloatVector.fromArray(SP, P8, i);
			final FloatVector p9 = FloatVector.fromArray(SP, P9, i);
			final FloatVector p10 = FloatVector.fromArray(SP, P10, i);
			final FloatVector p11 = FloatVector.fromArray(SP, P11, i);
			final FloatVector p12 = FloatVector.fromArray(SP, P12, i);
			final FloatVector t0 = p0.neg().add(s0).mul(c0);
			final FloatVector t1 = p1.neg().add(s1).mul(c1);
			final FloatVector t2 = p2.neg().add(s2).mul(c2);
			final FloatVector t3 = p3.neg().add(s3).mul(c3);
			final FloatVector t4 = p4.neg().add(s4).mul(c4);
			final FloatVector t5 = p5.neg().add(s5).mul(c5);
			final FloatVector t6 = p6.neg().add(s6).mul(c6);
			final FloatVector t7 = p7.neg().add(s7).mul(c7);
			final FloatVector t8 = p8.neg().add(s8).mul(c8);
			final FloatVector t9 = p9.neg().add(s9).mul(c9);
			final FloatVector t10 = p10.neg().add(s10).mul(c10);
			final FloatVector t11 = p11.neg().add(s11).mul(c11);
			final FloatVector t12 = p12.neg().add(s12).mul(c12);
			final FloatVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2)).add(p3.mul(t3)).add(p4.mul(t4)).add(p5.mul(t5)).add(p6.mul(t6)).add(p7.mul(t7)).add(p8.mul(t8)).add(p9.mul(t9)).add(p10.mul(t10)).add(p11.mul(t11)).add(p12.mul(t12));
			accCoal = accCoal.add(sc);
			final FloatVector d0 = p0.mul(sc.sub(t0)).add(p1.mul(rt[12])).add(p2.mul(rt[24])).add(p3.mul(rt[36])).add(p4.mul(rt[48])).add(p5.mul(rt[60])).add(p6.mul(rt[72])).add(p7.mul(rt[84])).add(p8.mul(rt[96])).add(p9.mul(rt[108])).add(p10.mul(rt[120])).add(p11.mul(rt[132])).add(p12.mul(rt[144])).sub(p0.mul(out[0]));
			d0.intoArray(D0, i); aD0 = aD0.add(d0);
			final FloatVector d1 = p1.mul(sc.sub(t1)).add(p0.mul(rt[0])).add(p2.mul(rt[25])).add(p3.mul(rt[37])).add(p4.mul(rt[49])).add(p5.mul(rt[61])).add(p6.mul(rt[73])).add(p7.mul(rt[85])).add(p8.mul(rt[97])).add(p9.mul(rt[109])).add(p10.mul(rt[121])).add(p11.mul(rt[133])).add(p12.mul(rt[145])).sub(p1.mul(out[1]));
			d1.intoArray(D1, i); aD1 = aD1.add(d1);
			final FloatVector d2 = p2.mul(sc.sub(t2)).add(p0.mul(rt[1])).add(p1.mul(rt[13])).add(p3.mul(rt[38])).add(p4.mul(rt[50])).add(p5.mul(rt[62])).add(p6.mul(rt[74])).add(p7.mul(rt[86])).add(p8.mul(rt[98])).add(p9.mul(rt[110])).add(p10.mul(rt[122])).add(p11.mul(rt[134])).add(p12.mul(rt[146])).sub(p2.mul(out[2]));
			d2.intoArray(D2, i); aD2 = aD2.add(d2);
			final FloatVector d3 = p3.mul(sc.sub(t3)).add(p0.mul(rt[2])).add(p1.mul(rt[14])).add(p2.mul(rt[26])).add(p4.mul(rt[51])).add(p5.mul(rt[63])).add(p6.mul(rt[75])).add(p7.mul(rt[87])).add(p8.mul(rt[99])).add(p9.mul(rt[111])).add(p10.mul(rt[123])).add(p11.mul(rt[135])).add(p12.mul(rt[147])).sub(p3.mul(out[3]));
			d3.intoArray(D3, i); aD3 = aD3.add(d3);
			final FloatVector d4 = p4.mul(sc.sub(t4)).add(p0.mul(rt[3])).add(p1.mul(rt[15])).add(p2.mul(rt[27])).add(p3.mul(rt[39])).add(p5.mul(rt[64])).add(p6.mul(rt[76])).add(p7.mul(rt[88])).add(p8.mul(rt[100])).add(p9.mul(rt[112])).add(p10.mul(rt[124])).add(p11.mul(rt[136])).add(p12.mul(rt[148])).sub(p4.mul(out[4]));
			d4.intoArray(D4, i); aD4 = aD4.add(d4);
			final FloatVector d5 = p5.mul(sc.sub(t5)).add(p0.mul(rt[4])).add(p1.mul(rt[16])).add(p2.mul(rt[28])).add(p3.mul(rt[40])).add(p4.mul(rt[52])).add(p6.mul(rt[77])).add(p7.mul(rt[89])).add(p8.mul(rt[101])).add(p9.mul(rt[113])).add(p10.mul(rt[125])).add(p11.mul(rt[137])).add(p12.mul(rt[149])).sub(p5.mul(out[5]));
			d5.intoArray(D5, i); aD5 = aD5.add(d5);
			final FloatVector d6 = p6.mul(sc.sub(t6)).add(p0.mul(rt[5])).add(p1.mul(rt[17])).add(p2.mul(rt[29])).add(p3.mul(rt[41])).add(p4.mul(rt[53])).add(p5.mul(rt[65])).add(p7.mul(rt[90])).add(p8.mul(rt[102])).add(p9.mul(rt[114])).add(p10.mul(rt[126])).add(p11.mul(rt[138])).add(p12.mul(rt[150])).sub(p6.mul(out[6]));
			d6.intoArray(D6, i); aD6 = aD6.add(d6);
			final FloatVector d7 = p7.mul(sc.sub(t7)).add(p0.mul(rt[6])).add(p1.mul(rt[18])).add(p2.mul(rt[30])).add(p3.mul(rt[42])).add(p4.mul(rt[54])).add(p5.mul(rt[66])).add(p6.mul(rt[78])).add(p8.mul(rt[103])).add(p9.mul(rt[115])).add(p10.mul(rt[127])).add(p11.mul(rt[139])).add(p12.mul(rt[151])).sub(p7.mul(out[7]));
			d7.intoArray(D7, i); aD7 = aD7.add(d7);
			final FloatVector d8 = p8.mul(sc.sub(t8)).add(p0.mul(rt[7])).add(p1.mul(rt[19])).add(p2.mul(rt[31])).add(p3.mul(rt[43])).add(p4.mul(rt[55])).add(p5.mul(rt[67])).add(p6.mul(rt[79])).add(p7.mul(rt[91])).add(p9.mul(rt[116])).add(p10.mul(rt[128])).add(p11.mul(rt[140])).add(p12.mul(rt[152])).sub(p8.mul(out[8]));
			d8.intoArray(D8, i); aD8 = aD8.add(d8);
			final FloatVector d9 = p9.mul(sc.sub(t9)).add(p0.mul(rt[8])).add(p1.mul(rt[20])).add(p2.mul(rt[32])).add(p3.mul(rt[44])).add(p4.mul(rt[56])).add(p5.mul(rt[68])).add(p6.mul(rt[80])).add(p7.mul(rt[92])).add(p8.mul(rt[104])).add(p10.mul(rt[129])).add(p11.mul(rt[141])).add(p12.mul(rt[153])).sub(p9.mul(out[9]));
			d9.intoArray(D9, i); aD9 = aD9.add(d9);
			final FloatVector d10 = p10.mul(sc.sub(t10)).add(p0.mul(rt[9])).add(p1.mul(rt[21])).add(p2.mul(rt[33])).add(p3.mul(rt[45])).add(p4.mul(rt[57])).add(p5.mul(rt[69])).add(p6.mul(rt[81])).add(p7.mul(rt[93])).add(p8.mul(rt[105])).add(p9.mul(rt[117])).add(p11.mul(rt[142])).add(p12.mul(rt[154])).sub(p10.mul(out[10]));
			d10.intoArray(D10, i); aD10 = aD10.add(d10);
			final FloatVector d11 = p11.mul(sc.sub(t11)).add(p0.mul(rt[10])).add(p1.mul(rt[22])).add(p2.mul(rt[34])).add(p3.mul(rt[46])).add(p4.mul(rt[58])).add(p5.mul(rt[70])).add(p6.mul(rt[82])).add(p7.mul(rt[94])).add(p8.mul(rt[106])).add(p9.mul(rt[118])).add(p10.mul(rt[130])).add(p12.mul(rt[155])).sub(p11.mul(out[11]));
			d11.intoArray(D11, i); aD11 = aD11.add(d11);
			final FloatVector d12 = p12.mul(sc.sub(t12)).add(p0.mul(rt[11])).add(p1.mul(rt[23])).add(p2.mul(rt[35])).add(p3.mul(rt[47])).add(p4.mul(rt[59])).add(p5.mul(rt[71])).add(p6.mul(rt[83])).add(p7.mul(rt[95])).add(p8.mul(rt[107])).add(p9.mul(rt[119])).add(p10.mul(rt[131])).add(p11.mul(rt[143])).sub(p12.mul(out[12]));
			d12.intoArray(D12, i); aD12 = aD12.add(d12);
		}
		dLast = -accCoal.reduceLanes(VectorOperators.ADD) / 2;
		sD[0] = aD0.reduceLanes(VectorOperators.ADD);
		sD[1] = aD1.reduceLanes(VectorOperators.ADD);
		sD[2] = aD2.reduceLanes(VectorOperators.ADD);
		sD[3] = aD3.reduceLanes(VectorOperators.ADD);
		sD[4] = aD4.reduceLanes(VectorOperators.ADD);
		sD[5] = aD5.reduceLanes(VectorOperators.ADD);
		sD[6] = aD6.reduceLanes(VectorOperators.ADD);
		sD[7] = aD7.reduceLanes(VectorOperators.ADD);
		sD[8] = aD8.reduceLanes(VectorOperators.ADD);
		sD[9] = aD9.reduceLanes(VectorOperators.ADD);
		sD[10] = aD10.reduceLanes(VectorOperators.ADD);
		sD[11] = aD11.reduceLanes(VectorOperators.ADD);
		sD[12] = aD12.reduceLanes(VectorOperators.ADD);
	}

	/** second derivative and the max of the approximate third derivative */
	private void pass2(final int Lpad) {
		final float[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12;
		final float[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4, D5 = this.D5, D6 = this.D6, D7 = this.D7, D8 = this.D8, D9 = this.D9, D10 = this.D10, D11 = this.D11, D12 = this.D12;
		final float[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4, DD5 = this.DD5, DD6 = this.DD6, DD7 = this.DD7, DD8 = this.DD8, DD9 = this.DD9, DD10 = this.DD10, DD11 = this.DD11, DD12 = this.DD12;
		final float[] rt = this.rt, out = this.out;
		final float s0 = (float) sS[0], s1 = (float) sS[1], s2 = (float) sS[2], s3 = (float) sS[3], s4 = (float) sS[4], s5 = (float) sS[5], s6 = (float) sS[6], s7 = (float) sS[7], s8 = (float) sS[8], s9 = (float) sS[9], s10 = (float) sS[10], s11 = (float) sS[11], s12 = (float) sS[12];
		final float c0 = cr[0], c1 = cr[1], c2 = cr[2], c3 = cr[3], c4 = cr[4], c5 = cr[5], c6 = cr[6], c7 = cr[7], c8 = cr[8], c9 = cr[9], c10 = cr[10], c11 = cr[11], c12 = cr[12];
		final float q0 = (float) sD[0], q1 = (float) sD[1], q2 = (float) sD[2], q3 = (float) sD[3], q4 = (float) sD[4], q5 = (float) sD[5], q6 = (float) sD[6], q7 = (float) sD[7], q8 = (float) sD[8], q9 = (float) sD[9], q10 = (float) sD[10], q11 = (float) sD[11], q12 = (float) sD[12];
		FloatVector accPc = FloatVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final FloatVector p0 = FloatVector.fromArray(SP, P0, i);
			final FloatVector p1 = FloatVector.fromArray(SP, P1, i);
			final FloatVector p2 = FloatVector.fromArray(SP, P2, i);
			final FloatVector p3 = FloatVector.fromArray(SP, P3, i);
			final FloatVector p4 = FloatVector.fromArray(SP, P4, i);
			final FloatVector p5 = FloatVector.fromArray(SP, P5, i);
			final FloatVector p6 = FloatVector.fromArray(SP, P6, i);
			final FloatVector p7 = FloatVector.fromArray(SP, P7, i);
			final FloatVector p8 = FloatVector.fromArray(SP, P8, i);
			final FloatVector p9 = FloatVector.fromArray(SP, P9, i);
			final FloatVector p10 = FloatVector.fromArray(SP, P10, i);
			final FloatVector p11 = FloatVector.fromArray(SP, P11, i);
			final FloatVector p12 = FloatVector.fromArray(SP, P12, i);
			final FloatVector d0 = FloatVector.fromArray(SP, D0, i);
			final FloatVector d1 = FloatVector.fromArray(SP, D1, i);
			final FloatVector d2 = FloatVector.fromArray(SP, D2, i);
			final FloatVector d3 = FloatVector.fromArray(SP, D3, i);
			final FloatVector d4 = FloatVector.fromArray(SP, D4, i);
			final FloatVector d5 = FloatVector.fromArray(SP, D5, i);
			final FloatVector d6 = FloatVector.fromArray(SP, D6, i);
			final FloatVector d7 = FloatVector.fromArray(SP, D7, i);
			final FloatVector d8 = FloatVector.fromArray(SP, D8, i);
			final FloatVector d9 = FloatVector.fromArray(SP, D9, i);
			final FloatVector d10 = FloatVector.fromArray(SP, D10, i);
			final FloatVector d11 = FloatVector.fromArray(SP, D11, i);
			final FloatVector d12 = FloatVector.fromArray(SP, D12, i);
			final FloatVector t0 = p0.neg().add(s0).mul(c0);
			final FloatVector t1 = p1.neg().add(s1).mul(c1);
			final FloatVector t2 = p2.neg().add(s2).mul(c2);
			final FloatVector t3 = p3.neg().add(s3).mul(c3);
			final FloatVector t4 = p4.neg().add(s4).mul(c4);
			final FloatVector t5 = p5.neg().add(s5).mul(c5);
			final FloatVector t6 = p6.neg().add(s6).mul(c6);
			final FloatVector t7 = p7.neg().add(s7).mul(c7);
			final FloatVector t8 = p8.neg().add(s8).mul(c8);
			final FloatVector t9 = p9.neg().add(s9).mul(c9);
			final FloatVector t10 = p10.neg().add(s10).mul(c10);
			final FloatVector t11 = p11.neg().add(s11).mul(c11);
			final FloatVector t12 = p12.neg().add(s12).mul(c12);
			final FloatVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2)).add(p3.mul(t3)).add(p4.mul(t4)).add(p5.mul(t5)).add(p6.mul(t6)).add(p7.mul(t7)).add(p8.mul(t8)).add(p9.mul(t9)).add(p10.mul(t10)).add(p11.mul(t11)).add(p12.mul(t12));
			final FloatVector pc = d0.mul(p0.mul(-2.0f).add(s0)).add(p0.mul(q0)).mul(c0).add(d1.mul(p1.mul(-2.0f).add(s1)).add(p1.mul(q1)).mul(c1)).add(d2.mul(p2.mul(-2.0f).add(s2)).add(p2.mul(q2)).mul(c2)).add(d3.mul(p3.mul(-2.0f).add(s3)).add(p3.mul(q3)).mul(c3)).add(d4.mul(p4.mul(-2.0f).add(s4)).add(p4.mul(q4)).mul(c4)).add(d5.mul(p5.mul(-2.0f).add(s5)).add(p5.mul(q5)).mul(c5)).add(d6.mul(p6.mul(-2.0f).add(s6)).add(p6.mul(q6)).mul(c6)).add(d7.mul(p7.mul(-2.0f).add(s7)).add(p7.mul(q7)).mul(c7)).add(d8.mul(p8.mul(-2.0f).add(s8)).add(p8.mul(q8)).mul(c8)).add(d9.mul(p9.mul(-2.0f).add(s9)).add(p9.mul(q9)).mul(c9)).add(d10.mul(p10.mul(-2.0f).add(s10)).add(p10.mul(q10)).mul(c10)).add(d11.mul(p11.mul(-2.0f).add(s11)).add(p11.mul(q11)).mul(c11)).add(d12.mul(p12.mul(-2.0f).add(s12)).add(p12.mul(q12)).mul(c12));
			accPc = accPc.add(pc);
			final FloatVector dd0 = sc.sub(t0).mul(d0).add(p0.mul(pc.sub(d0.neg().add(q0).mul(c0)))).add(d1.mul(rt[12])).add(d2.mul(rt[24])).add(d3.mul(rt[36])).add(d4.mul(rt[48])).add(d5.mul(rt[60])).add(d6.mul(rt[72])).add(d7.mul(rt[84])).add(d8.mul(rt[96])).add(d9.mul(rt[108])).add(d10.mul(rt[120])).add(d11.mul(rt[132])).add(d12.mul(rt[144])).sub(d0.mul(out[0]));
			dd0.intoArray(DD0, i);
			final FloatVector dd1 = sc.sub(t1).mul(d1).add(p1.mul(pc.sub(d1.neg().add(q1).mul(c1)))).add(d0.mul(rt[0])).add(d2.mul(rt[25])).add(d3.mul(rt[37])).add(d4.mul(rt[49])).add(d5.mul(rt[61])).add(d6.mul(rt[73])).add(d7.mul(rt[85])).add(d8.mul(rt[97])).add(d9.mul(rt[109])).add(d10.mul(rt[121])).add(d11.mul(rt[133])).add(d12.mul(rt[145])).sub(d1.mul(out[1]));
			dd1.intoArray(DD1, i);
			final FloatVector dd2 = sc.sub(t2).mul(d2).add(p2.mul(pc.sub(d2.neg().add(q2).mul(c2)))).add(d0.mul(rt[1])).add(d1.mul(rt[13])).add(d3.mul(rt[38])).add(d4.mul(rt[50])).add(d5.mul(rt[62])).add(d6.mul(rt[74])).add(d7.mul(rt[86])).add(d8.mul(rt[98])).add(d9.mul(rt[110])).add(d10.mul(rt[122])).add(d11.mul(rt[134])).add(d12.mul(rt[146])).sub(d2.mul(out[2]));
			dd2.intoArray(DD2, i);
			final FloatVector dd3 = sc.sub(t3).mul(d3).add(p3.mul(pc.sub(d3.neg().add(q3).mul(c3)))).add(d0.mul(rt[2])).add(d1.mul(rt[14])).add(d2.mul(rt[26])).add(d4.mul(rt[51])).add(d5.mul(rt[63])).add(d6.mul(rt[75])).add(d7.mul(rt[87])).add(d8.mul(rt[99])).add(d9.mul(rt[111])).add(d10.mul(rt[123])).add(d11.mul(rt[135])).add(d12.mul(rt[147])).sub(d3.mul(out[3]));
			dd3.intoArray(DD3, i);
			final FloatVector dd4 = sc.sub(t4).mul(d4).add(p4.mul(pc.sub(d4.neg().add(q4).mul(c4)))).add(d0.mul(rt[3])).add(d1.mul(rt[15])).add(d2.mul(rt[27])).add(d3.mul(rt[39])).add(d5.mul(rt[64])).add(d6.mul(rt[76])).add(d7.mul(rt[88])).add(d8.mul(rt[100])).add(d9.mul(rt[112])).add(d10.mul(rt[124])).add(d11.mul(rt[136])).add(d12.mul(rt[148])).sub(d4.mul(out[4]));
			dd4.intoArray(DD4, i);
			final FloatVector dd5 = sc.sub(t5).mul(d5).add(p5.mul(pc.sub(d5.neg().add(q5).mul(c5)))).add(d0.mul(rt[4])).add(d1.mul(rt[16])).add(d2.mul(rt[28])).add(d3.mul(rt[40])).add(d4.mul(rt[52])).add(d6.mul(rt[77])).add(d7.mul(rt[89])).add(d8.mul(rt[101])).add(d9.mul(rt[113])).add(d10.mul(rt[125])).add(d11.mul(rt[137])).add(d12.mul(rt[149])).sub(d5.mul(out[5]));
			dd5.intoArray(DD5, i);
			final FloatVector dd6 = sc.sub(t6).mul(d6).add(p6.mul(pc.sub(d6.neg().add(q6).mul(c6)))).add(d0.mul(rt[5])).add(d1.mul(rt[17])).add(d2.mul(rt[29])).add(d3.mul(rt[41])).add(d4.mul(rt[53])).add(d5.mul(rt[65])).add(d7.mul(rt[90])).add(d8.mul(rt[102])).add(d9.mul(rt[114])).add(d10.mul(rt[126])).add(d11.mul(rt[138])).add(d12.mul(rt[150])).sub(d6.mul(out[6]));
			dd6.intoArray(DD6, i);
			final FloatVector dd7 = sc.sub(t7).mul(d7).add(p7.mul(pc.sub(d7.neg().add(q7).mul(c7)))).add(d0.mul(rt[6])).add(d1.mul(rt[18])).add(d2.mul(rt[30])).add(d3.mul(rt[42])).add(d4.mul(rt[54])).add(d5.mul(rt[66])).add(d6.mul(rt[78])).add(d8.mul(rt[103])).add(d9.mul(rt[115])).add(d10.mul(rt[127])).add(d11.mul(rt[139])).add(d12.mul(rt[151])).sub(d7.mul(out[7]));
			dd7.intoArray(DD7, i);
			final FloatVector dd8 = sc.sub(t8).mul(d8).add(p8.mul(pc.sub(d8.neg().add(q8).mul(c8)))).add(d0.mul(rt[7])).add(d1.mul(rt[19])).add(d2.mul(rt[31])).add(d3.mul(rt[43])).add(d4.mul(rt[55])).add(d5.mul(rt[67])).add(d6.mul(rt[79])).add(d7.mul(rt[91])).add(d9.mul(rt[116])).add(d10.mul(rt[128])).add(d11.mul(rt[140])).add(d12.mul(rt[152])).sub(d8.mul(out[8]));
			dd8.intoArray(DD8, i);
			final FloatVector dd9 = sc.sub(t9).mul(d9).add(p9.mul(pc.sub(d9.neg().add(q9).mul(c9)))).add(d0.mul(rt[8])).add(d1.mul(rt[20])).add(d2.mul(rt[32])).add(d3.mul(rt[44])).add(d4.mul(rt[56])).add(d5.mul(rt[68])).add(d6.mul(rt[80])).add(d7.mul(rt[92])).add(d8.mul(rt[104])).add(d10.mul(rt[129])).add(d11.mul(rt[141])).add(d12.mul(rt[153])).sub(d9.mul(out[9]));
			dd9.intoArray(DD9, i);
			final FloatVector dd10 = sc.sub(t10).mul(d10).add(p10.mul(pc.sub(d10.neg().add(q10).mul(c10)))).add(d0.mul(rt[9])).add(d1.mul(rt[21])).add(d2.mul(rt[33])).add(d3.mul(rt[45])).add(d4.mul(rt[57])).add(d5.mul(rt[69])).add(d6.mul(rt[81])).add(d7.mul(rt[93])).add(d8.mul(rt[105])).add(d9.mul(rt[117])).add(d11.mul(rt[142])).add(d12.mul(rt[154])).sub(d10.mul(out[10]));
			dd10.intoArray(DD10, i);
			final FloatVector dd11 = sc.sub(t11).mul(d11).add(p11.mul(pc.sub(d11.neg().add(q11).mul(c11)))).add(d0.mul(rt[10])).add(d1.mul(rt[22])).add(d2.mul(rt[34])).add(d3.mul(rt[46])).add(d4.mul(rt[58])).add(d5.mul(rt[70])).add(d6.mul(rt[82])).add(d7.mul(rt[94])).add(d8.mul(rt[106])).add(d9.mul(rt[118])).add(d10.mul(rt[130])).add(d12.mul(rt[155])).sub(d11.mul(out[11]));
			dd11.intoArray(DD11, i);
			final FloatVector dd12 = sc.sub(t12).mul(d12).add(p12.mul(pc.sub(d12.neg().add(q12).mul(c12)))).add(d0.mul(rt[11])).add(d1.mul(rt[23])).add(d2.mul(rt[35])).add(d3.mul(rt[47])).add(d4.mul(rt[59])).add(d5.mul(rt[71])).add(d6.mul(rt[83])).add(d7.mul(rt[95])).add(d8.mul(rt[107])).add(d9.mul(rt[119])).add(d10.mul(rt[131])).add(d11.mul(rt[143])).sub(d12.mul(out[12]));
			dd12.intoArray(DD12, i);
		}
		ddLast = -accPc.reduceLanes(VectorOperators.ADD) / 2;
	}

	/** max of the approximate third derivative */
	private void pass2b(final int Lpad) {
		final float[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12;
		final float[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4, DD5 = this.DD5, DD6 = this.DD6, DD7 = this.DD7, DD8 = this.DD8, DD9 = this.DD9, DD10 = this.DD10, DD11 = this.DD11, DD12 = this.DD12;
		final float[] rt = this.rt, out = this.out;
		final float s0 = (float) sS[0], s1 = (float) sS[1], s2 = (float) sS[2], s3 = (float) sS[3], s4 = (float) sS[4], s5 = (float) sS[5], s6 = (float) sS[6], s7 = (float) sS[7], s8 = (float) sS[8], s9 = (float) sS[9], s10 = (float) sS[10], s11 = (float) sS[11], s12 = (float) sS[12];
		final float c0 = cr[0], c1 = cr[1], c2 = cr[2], c3 = cr[3], c4 = cr[4], c5 = cr[5], c6 = cr[6], c7 = cr[7], c8 = cr[8], c9 = cr[9], c10 = cr[10], c11 = cr[11], c12 = cr[12];
		FloatVector accMax = FloatVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final FloatVector p0 = FloatVector.fromArray(SP, P0, i);
			final FloatVector p1 = FloatVector.fromArray(SP, P1, i);
			final FloatVector p2 = FloatVector.fromArray(SP, P2, i);
			final FloatVector p3 = FloatVector.fromArray(SP, P3, i);
			final FloatVector p4 = FloatVector.fromArray(SP, P4, i);
			final FloatVector p5 = FloatVector.fromArray(SP, P5, i);
			final FloatVector p6 = FloatVector.fromArray(SP, P6, i);
			final FloatVector p7 = FloatVector.fromArray(SP, P7, i);
			final FloatVector p8 = FloatVector.fromArray(SP, P8, i);
			final FloatVector p9 = FloatVector.fromArray(SP, P9, i);
			final FloatVector p10 = FloatVector.fromArray(SP, P10, i);
			final FloatVector p11 = FloatVector.fromArray(SP, P11, i);
			final FloatVector p12 = FloatVector.fromArray(SP, P12, i);
			final FloatVector dd0 = FloatVector.fromArray(SP, DD0, i);
			final FloatVector dd1 = FloatVector.fromArray(SP, DD1, i);
			final FloatVector dd2 = FloatVector.fromArray(SP, DD2, i);
			final FloatVector dd3 = FloatVector.fromArray(SP, DD3, i);
			final FloatVector dd4 = FloatVector.fromArray(SP, DD4, i);
			final FloatVector dd5 = FloatVector.fromArray(SP, DD5, i);
			final FloatVector dd6 = FloatVector.fromArray(SP, DD6, i);
			final FloatVector dd7 = FloatVector.fromArray(SP, DD7, i);
			final FloatVector dd8 = FloatVector.fromArray(SP, DD8, i);
			final FloatVector dd9 = FloatVector.fromArray(SP, DD9, i);
			final FloatVector dd10 = FloatVector.fromArray(SP, DD10, i);
			final FloatVector dd11 = FloatVector.fromArray(SP, DD11, i);
			final FloatVector dd12 = FloatVector.fromArray(SP, DD12, i);
			final FloatVector t0 = p0.neg().add(s0).mul(c0);
			final FloatVector t1 = p1.neg().add(s1).mul(c1);
			final FloatVector t2 = p2.neg().add(s2).mul(c2);
			final FloatVector t3 = p3.neg().add(s3).mul(c3);
			final FloatVector t4 = p4.neg().add(s4).mul(c4);
			final FloatVector t5 = p5.neg().add(s5).mul(c5);
			final FloatVector t6 = p6.neg().add(s6).mul(c6);
			final FloatVector t7 = p7.neg().add(s7).mul(c7);
			final FloatVector t8 = p8.neg().add(s8).mul(c8);
			final FloatVector t9 = p9.neg().add(s9).mul(c9);
			final FloatVector t10 = p10.neg().add(s10).mul(c10);
			final FloatVector t11 = p11.neg().add(s11).mul(c11);
			final FloatVector t12 = p12.neg().add(s12).mul(c12);
			final FloatVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2)).add(p3.mul(t3)).add(p4.mul(t4)).add(p5.mul(t5)).add(p6.mul(t6)).add(p7.mul(t7)).add(p8.mul(t8)).add(p9.mul(t9)).add(p10.mul(t10)).add(p11.mul(t11)).add(p12.mul(t12));
			accMax = accMax.max(sc.sub(t0).mul(dd0).add(dd1.mul(rt[12])).add(dd2.mul(rt[24])).add(dd3.mul(rt[36])).add(dd4.mul(rt[48])).add(dd5.mul(rt[60])).add(dd6.mul(rt[72])).add(dd7.mul(rt[84])).add(dd8.mul(rt[96])).add(dd9.mul(rt[108])).add(dd10.mul(rt[120])).add(dd11.mul(rt[132])).add(dd12.mul(rt[144])).sub(dd0.mul(out[0])).abs());
			accMax = accMax.max(sc.sub(t1).mul(dd1).add(dd0.mul(rt[0])).add(dd2.mul(rt[25])).add(dd3.mul(rt[37])).add(dd4.mul(rt[49])).add(dd5.mul(rt[61])).add(dd6.mul(rt[73])).add(dd7.mul(rt[85])).add(dd8.mul(rt[97])).add(dd9.mul(rt[109])).add(dd10.mul(rt[121])).add(dd11.mul(rt[133])).add(dd12.mul(rt[145])).sub(dd1.mul(out[1])).abs());
			accMax = accMax.max(sc.sub(t2).mul(dd2).add(dd0.mul(rt[1])).add(dd1.mul(rt[13])).add(dd3.mul(rt[38])).add(dd4.mul(rt[50])).add(dd5.mul(rt[62])).add(dd6.mul(rt[74])).add(dd7.mul(rt[86])).add(dd8.mul(rt[98])).add(dd9.mul(rt[110])).add(dd10.mul(rt[122])).add(dd11.mul(rt[134])).add(dd12.mul(rt[146])).sub(dd2.mul(out[2])).abs());
			accMax = accMax.max(sc.sub(t3).mul(dd3).add(dd0.mul(rt[2])).add(dd1.mul(rt[14])).add(dd2.mul(rt[26])).add(dd4.mul(rt[51])).add(dd5.mul(rt[63])).add(dd6.mul(rt[75])).add(dd7.mul(rt[87])).add(dd8.mul(rt[99])).add(dd9.mul(rt[111])).add(dd10.mul(rt[123])).add(dd11.mul(rt[135])).add(dd12.mul(rt[147])).sub(dd3.mul(out[3])).abs());
			accMax = accMax.max(sc.sub(t4).mul(dd4).add(dd0.mul(rt[3])).add(dd1.mul(rt[15])).add(dd2.mul(rt[27])).add(dd3.mul(rt[39])).add(dd5.mul(rt[64])).add(dd6.mul(rt[76])).add(dd7.mul(rt[88])).add(dd8.mul(rt[100])).add(dd9.mul(rt[112])).add(dd10.mul(rt[124])).add(dd11.mul(rt[136])).add(dd12.mul(rt[148])).sub(dd4.mul(out[4])).abs());
			accMax = accMax.max(sc.sub(t5).mul(dd5).add(dd0.mul(rt[4])).add(dd1.mul(rt[16])).add(dd2.mul(rt[28])).add(dd3.mul(rt[40])).add(dd4.mul(rt[52])).add(dd6.mul(rt[77])).add(dd7.mul(rt[89])).add(dd8.mul(rt[101])).add(dd9.mul(rt[113])).add(dd10.mul(rt[125])).add(dd11.mul(rt[137])).add(dd12.mul(rt[149])).sub(dd5.mul(out[5])).abs());
			accMax = accMax.max(sc.sub(t6).mul(dd6).add(dd0.mul(rt[5])).add(dd1.mul(rt[17])).add(dd2.mul(rt[29])).add(dd3.mul(rt[41])).add(dd4.mul(rt[53])).add(dd5.mul(rt[65])).add(dd7.mul(rt[90])).add(dd8.mul(rt[102])).add(dd9.mul(rt[114])).add(dd10.mul(rt[126])).add(dd11.mul(rt[138])).add(dd12.mul(rt[150])).sub(dd6.mul(out[6])).abs());
			accMax = accMax.max(sc.sub(t7).mul(dd7).add(dd0.mul(rt[6])).add(dd1.mul(rt[18])).add(dd2.mul(rt[30])).add(dd3.mul(rt[42])).add(dd4.mul(rt[54])).add(dd5.mul(rt[66])).add(dd6.mul(rt[78])).add(dd8.mul(rt[103])).add(dd9.mul(rt[115])).add(dd10.mul(rt[127])).add(dd11.mul(rt[139])).add(dd12.mul(rt[151])).sub(dd7.mul(out[7])).abs());
			accMax = accMax.max(sc.sub(t8).mul(dd8).add(dd0.mul(rt[7])).add(dd1.mul(rt[19])).add(dd2.mul(rt[31])).add(dd3.mul(rt[43])).add(dd4.mul(rt[55])).add(dd5.mul(rt[67])).add(dd6.mul(rt[79])).add(dd7.mul(rt[91])).add(dd9.mul(rt[116])).add(dd10.mul(rt[128])).add(dd11.mul(rt[140])).add(dd12.mul(rt[152])).sub(dd8.mul(out[8])).abs());
			accMax = accMax.max(sc.sub(t9).mul(dd9).add(dd0.mul(rt[8])).add(dd1.mul(rt[20])).add(dd2.mul(rt[32])).add(dd3.mul(rt[44])).add(dd4.mul(rt[56])).add(dd5.mul(rt[68])).add(dd6.mul(rt[80])).add(dd7.mul(rt[92])).add(dd8.mul(rt[104])).add(dd10.mul(rt[129])).add(dd11.mul(rt[141])).add(dd12.mul(rt[153])).sub(dd9.mul(out[9])).abs());
			accMax = accMax.max(sc.sub(t10).mul(dd10).add(dd0.mul(rt[9])).add(dd1.mul(rt[21])).add(dd2.mul(rt[33])).add(dd3.mul(rt[45])).add(dd4.mul(rt[57])).add(dd5.mul(rt[69])).add(dd6.mul(rt[81])).add(dd7.mul(rt[93])).add(dd8.mul(rt[105])).add(dd9.mul(rt[117])).add(dd11.mul(rt[142])).add(dd12.mul(rt[154])).sub(dd10.mul(out[10])).abs());
			accMax = accMax.max(sc.sub(t11).mul(dd11).add(dd0.mul(rt[10])).add(dd1.mul(rt[22])).add(dd2.mul(rt[34])).add(dd3.mul(rt[46])).add(dd4.mul(rt[58])).add(dd5.mul(rt[70])).add(dd6.mul(rt[82])).add(dd7.mul(rt[94])).add(dd8.mul(rt[106])).add(dd9.mul(rt[118])).add(dd10.mul(rt[130])).add(dd12.mul(rt[155])).sub(dd11.mul(out[11])).abs());
			accMax = accMax.max(sc.sub(t12).mul(dd12).add(dd0.mul(rt[11])).add(dd1.mul(rt[23])).add(dd2.mul(rt[35])).add(dd3.mul(rt[47])).add(dd4.mul(rt[59])).add(dd5.mul(rt[71])).add(dd6.mul(rt[83])).add(dd7.mul(rt[95])).add(dd8.mul(rt[107])).add(dd9.mul(rt[119])).add(dd10.mul(rt[131])).add(dd11.mul(rt[143])).sub(dd12.mul(out[12])).abs());
		}
		maxDDD = accMax.reduceLanes(VectorOperators.MAX);
	}

	/** true if any probability would leave [0,1] or change by more than 0.2 with this step */
	private boolean needsSearch(final int Lpad, final float ts, final float tss) {
		final float[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12;
		final float[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4, D5 = this.D5, D6 = this.D6, D7 = this.D7, D8 = this.D8, D9 = this.D9, D10 = this.D10, D11 = this.D11, D12 = this.D12;
		final float[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4, DD5 = this.DD5, DD6 = this.DD6, DD7 = this.DD7, DD8 = this.DD8, DD9 = this.DD9, DD10 = this.DD10, DD11 = this.DD11, DD12 = this.DD12;
		for (int i = 0; i < Lpad; i += W) {
			VectorMask<Float> m = SP.maskAll(false);
			{
				final FloatVector p = FloatVector.fromArray(SP, P0, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D0, i).mul(ts)).add(FloatVector.fromArray(SP, DD0, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P1, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D1, i).mul(ts)).add(FloatVector.fromArray(SP, DD1, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P2, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D2, i).mul(ts)).add(FloatVector.fromArray(SP, DD2, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P3, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D3, i).mul(ts)).add(FloatVector.fromArray(SP, DD3, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P4, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D4, i).mul(ts)).add(FloatVector.fromArray(SP, DD4, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P5, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D5, i).mul(ts)).add(FloatVector.fromArray(SP, DD5, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P6, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D6, i).mul(ts)).add(FloatVector.fromArray(SP, DD6, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P7, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D7, i).mul(ts)).add(FloatVector.fromArray(SP, DD7, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P8, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D8, i).mul(ts)).add(FloatVector.fromArray(SP, DD8, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P9, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D9, i).mul(ts)).add(FloatVector.fromArray(SP, DD9, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P10, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D10, i).mul(ts)).add(FloatVector.fromArray(SP, DD10, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P11, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D11, i).mul(ts)).add(FloatVector.fromArray(SP, DD11, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			{
				final FloatVector p = FloatVector.fromArray(SP, P12, i);
				final FloatVector nv = p.add(FloatVector.fromArray(SP, D12, i).mul(ts)).add(FloatVector.fromArray(SP, DD12, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0f)).or(nv.compare(VectorOperators.LT, 0.0f)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2f));
			}
			if (m.anyTrue()) return true;
		}
		return false;
	}

	/** update and normalise; returns true if a probability became negative */
	private boolean pass4(final int Lpad, final int L, final float ts, final float tss) {
		final float[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12;
		final float[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4, D5 = this.D5, D6 = this.D6, D7 = this.D7, D8 = this.D8, D9 = this.D9, D10 = this.D10, D11 = this.D11, D12 = this.D12;
		final float[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4, DD5 = this.DD5, DD6 = this.DD6, DD7 = this.DD7, DD8 = this.DD8, DD9 = this.DD9, DD10 = this.DD10, DD11 = this.DD11, DD12 = this.DD12;
		final int lastBlock = Lpad - W;
		final VectorMask<Float> pad = SP.indexInRange(lastBlock, L).not();
		boolean failed = false;
		FloatVector aS0 = FloatVector.zero(SP), aS1 = FloatVector.zero(SP), aS2 = FloatVector.zero(SP), aS3 = FloatVector.zero(SP), aS4 = FloatVector.zero(SP), aS5 = FloatVector.zero(SP), aS6 = FloatVector.zero(SP), aS7 = FloatVector.zero(SP), aS8 = FloatVector.zero(SP), aS9 = FloatVector.zero(SP), aS10 = FloatVector.zero(SP), aS11 = FloatVector.zero(SP), aS12 = FloatVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final FloatVector x0 = FloatVector.fromArray(SP, P0, i).add(FloatVector.fromArray(SP, D0, i).mul(ts).add(FloatVector.fromArray(SP, DD0, i).mul(tss)));
			final FloatVector x1 = FloatVector.fromArray(SP, P1, i).add(FloatVector.fromArray(SP, D1, i).mul(ts).add(FloatVector.fromArray(SP, DD1, i).mul(tss)));
			final FloatVector x2 = FloatVector.fromArray(SP, P2, i).add(FloatVector.fromArray(SP, D2, i).mul(ts).add(FloatVector.fromArray(SP, DD2, i).mul(tss)));
			final FloatVector x3 = FloatVector.fromArray(SP, P3, i).add(FloatVector.fromArray(SP, D3, i).mul(ts).add(FloatVector.fromArray(SP, DD3, i).mul(tss)));
			final FloatVector x4 = FloatVector.fromArray(SP, P4, i).add(FloatVector.fromArray(SP, D4, i).mul(ts).add(FloatVector.fromArray(SP, DD4, i).mul(tss)));
			final FloatVector x5 = FloatVector.fromArray(SP, P5, i).add(FloatVector.fromArray(SP, D5, i).mul(ts).add(FloatVector.fromArray(SP, DD5, i).mul(tss)));
			final FloatVector x6 = FloatVector.fromArray(SP, P6, i).add(FloatVector.fromArray(SP, D6, i).mul(ts).add(FloatVector.fromArray(SP, DD6, i).mul(tss)));
			final FloatVector x7 = FloatVector.fromArray(SP, P7, i).add(FloatVector.fromArray(SP, D7, i).mul(ts).add(FloatVector.fromArray(SP, DD7, i).mul(tss)));
			final FloatVector x8 = FloatVector.fromArray(SP, P8, i).add(FloatVector.fromArray(SP, D8, i).mul(ts).add(FloatVector.fromArray(SP, DD8, i).mul(tss)));
			final FloatVector x9 = FloatVector.fromArray(SP, P9, i).add(FloatVector.fromArray(SP, D9, i).mul(ts).add(FloatVector.fromArray(SP, DD9, i).mul(tss)));
			final FloatVector x10 = FloatVector.fromArray(SP, P10, i).add(FloatVector.fromArray(SP, D10, i).mul(ts).add(FloatVector.fromArray(SP, DD10, i).mul(tss)));
			final FloatVector x11 = FloatVector.fromArray(SP, P11, i).add(FloatVector.fromArray(SP, D11, i).mul(ts).add(FloatVector.fromArray(SP, DD11, i).mul(tss)));
			final FloatVector x12 = FloatVector.fromArray(SP, P12, i).add(FloatVector.fromArray(SP, D12, i).mul(ts).add(FloatVector.fromArray(SP, DD12, i).mul(tss)));
			final VectorMask<Float> negative = x0.compare(VectorOperators.LT, 0.0f).or(x1.compare(VectorOperators.LT, 0.0f)).or(x2.compare(VectorOperators.LT, 0.0f)).or(x3.compare(VectorOperators.LT, 0.0f)).or(x4.compare(VectorOperators.LT, 0.0f)).or(x5.compare(VectorOperators.LT, 0.0f)).or(x6.compare(VectorOperators.LT, 0.0f)).or(x7.compare(VectorOperators.LT, 0.0f)).or(x8.compare(VectorOperators.LT, 0.0f)).or(x9.compare(VectorOperators.LT, 0.0f)).or(x10.compare(VectorOperators.LT, 0.0f)).or(x11.compare(VectorOperators.LT, 0.0f)).or(x12.compare(VectorOperators.LT, 0.0f));
			final FloatVector sum = x0.add(x1).add(x2).add(x3).add(x4).add(x5).add(x6).add(x7).add(x8).add(x9).add(x10).add(x11).add(x12);
			final boolean anyNegative = negative.anyTrue();
			if (anyNegative) failed = true;
			FloatVector y0 = x0.div(sum);
			FloatVector y1 = x1.div(sum);
			FloatVector y2 = x2.div(sum);
			FloatVector y3 = x3.div(sum);
			FloatVector y4 = x4.div(sum);
			FloatVector y5 = x5.div(sum);
			FloatVector y6 = x6.div(sum);
			FloatVector y7 = x7.div(sum);
			FloatVector y8 = x8.div(sum);
			FloatVector y9 = x9.div(sum);
			FloatVector y10 = x10.div(sum);
			FloatVector y11 = x11.div(sum);
			FloatVector y12 = x12.div(sum);
			if (anyNegative) {
				// like normalise(): leave the lineage unnormalised
				y0 = y0.blend(x0, negative);
				y1 = y1.blend(x1, negative);
				y2 = y2.blend(x2, negative);
				y3 = y3.blend(x3, negative);
				y4 = y4.blend(x4, negative);
				y5 = y5.blend(x5, negative);
				y6 = y6.blend(x6, negative);
				y7 = y7.blend(x7, negative);
				y8 = y8.blend(x8, negative);
				y9 = y9.blend(x9, negative);
				y10 = y10.blend(x10, negative);
				y11 = y11.blend(x11, negative);
				y12 = y12.blend(x12, negative);
			}
			if (i == lastBlock) {
				// padding lineages stay at zero (0/0 would give NaN)
				y0 = y0.blend(0.0f, pad);
				y1 = y1.blend(0.0f, pad);
				y2 = y2.blend(0.0f, pad);
				y3 = y3.blend(0.0f, pad);
				y4 = y4.blend(0.0f, pad);
				y5 = y5.blend(0.0f, pad);
				y6 = y6.blend(0.0f, pad);
				y7 = y7.blend(0.0f, pad);
				y8 = y8.blend(0.0f, pad);
				y9 = y9.blend(0.0f, pad);
				y10 = y10.blend(0.0f, pad);
				y11 = y11.blend(0.0f, pad);
				y12 = y12.blend(0.0f, pad);
			}
			y0.intoArray(P0, i); aS0 = aS0.add(y0);
			y1.intoArray(P1, i); aS1 = aS1.add(y1);
			y2.intoArray(P2, i); aS2 = aS2.add(y2);
			y3.intoArray(P3, i); aS3 = aS3.add(y3);
			y4.intoArray(P4, i); aS4 = aS4.add(y4);
			y5.intoArray(P5, i); aS5 = aS5.add(y5);
			y6.intoArray(P6, i); aS6 = aS6.add(y6);
			y7.intoArray(P7, i); aS7 = aS7.add(y7);
			y8.intoArray(P8, i); aS8 = aS8.add(y8);
			y9.intoArray(P9, i); aS9 = aS9.add(y9);
			y10.intoArray(P10, i); aS10 = aS10.add(y10);
			y11.intoArray(P11, i); aS11 = aS11.add(y11);
			y12.intoArray(P12, i); aS12 = aS12.add(y12);
		}
		sS[0] = aS0.reduceLanes(VectorOperators.ADD);
		sS[1] = aS1.reduceLanes(VectorOperators.ADD);
		sS[2] = aS2.reduceLanes(VectorOperators.ADD);
		sS[3] = aS3.reduceLanes(VectorOperators.ADD);
		sS[4] = aS4.reduceLanes(VectorOperators.ADD);
		sS[5] = aS5.reduceLanes(VectorOperators.ADD);
		sS[6] = aS6.reduceLanes(VectorOperators.ADD);
		sS[7] = aS7.reduceLanes(VectorOperators.ADD);
		sS[8] = aS8.reduceLanes(VectorOperators.ADD);
		sS[9] = aS9.reduceLanes(VectorOperators.ADD);
		sS[10] = aS10.reduceLanes(VectorOperators.ADD);
		sS[11] = aS11.reduceLanes(VectorOperators.ADD);
		sS[12] = aS12.reduceLanes(VectorOperators.ADD);
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
		final float[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12;
		for (int i = 0, k = 0; i < L; i++, k += N) {
			P0[i] = (float) p[k + 0];
			P1[i] = (float) p[k + 1];
			P2[i] = (float) p[k + 2];
			P3[i] = (float) p[k + 3];
			P4[i] = (float) p[k + 4];
			P5[i] = (float) p[k + 5];
			P6[i] = (float) p[k + 6];
			P7[i] = (float) p[k + 7];
			P8[i] = (float) p[k + 8];
			P9[i] = (float) p[k + 9];
			P10[i] = (float) p[k + 10];
			P11[i] = (float) p[k + 11];
			P12[i] = (float) p[k + 12];
		}
		for (int i = L; i < Lpad; i++) {
			P0[i] = 0; P1[i] = 0; P2[i] = 0; P3[i] = 0; P4[i] = 0; P5[i] = 0; P6[i] = 0; P7[i] = 0; P8[i] = 0; P9[i] = 0; P10[i] = 0; P11[i] = 0; P12[i] = 0;
		}
		for (int j = 0; j < N; j++) cr[j] = (float) coalescent_rates[j];
		final int last = length - 1;
		double logPSlot = p[last];
		{ double s = 0; for (int i = 0; i < L; i++) s += P0[i]; sS[0] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P1[i]; sS[1] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P2[i]; sS[2] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P3[i]; sS[3] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P4[i]; sS[4] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P5[i]; sS[5] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P6[i]; sS[6] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P7[i]; sS[7] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P8[i]; sS[8] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P9[i]; sS[9] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P10[i]; sS[10] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P11[i]; sS[11] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P12[i]; sS[12] = s; }

		while (duration > 0) {
			iterations++;
			pass1(Lpad);
			pass2(Lpad);
			pass2b(Lpad);
			double timeStep = FastMath.min(FastMath.cbrt(epsilon * 6 / maxDDD), FastMath.min(duration, max_step));
			double timeStepSquare = timeStep * timeStep * 0.5;
			if (needsSearch(Lpad, (float) timeStep, (float) timeStepSquare)) {
				// same sequential search as Euler2ndOrder.updateP, in the original element order
				final float[][] PA = {P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12}, DA = {D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12}, DDA = {DD0, DD1, DD2, DD3, DD4, DD5, DD6, DD7, DD8, DD9, DD10, DD11, DD12};
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
								P12[L - 1] = Float.NEGATIVE_INFINITY;
								break search;
							}
						}
					}
				}
			}
			if (P12[L - 1] == Float.NEGATIVE_INFINITY) {
				duration = 0.0;
			} else {
				logPSlot += dLast * timeStep + ddLast * timeStepSquare;
				if (pass4(Lpad, L, (float) timeStep, (float) timeStepSquare)) {
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
			p[k + 5] = P5[i];
			p[k + 6] = P6[i];
			p[k + 7] = P7[i];
			p[k + 8] = P8[i];
			p[k + 9] = P9[i];
			p[k + 10] = P10[i];
			p[k + 11] = P11[i];
			p[k + 12] = P12[i];
		}
		p[last] = logPSlot;
	}
}
