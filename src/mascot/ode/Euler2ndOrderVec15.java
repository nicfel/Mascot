package mascot.ode;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.apache.commons.math3.util.FastMath;

/**
 * SIMD (Vector API) solver for 15 states in double precision, written out like Euler2ndOrderVec3:
 * several lineages are processed per instruction, per-state values stay in vector registers and
 * all migration terms are unrolled. Lineage state probabilities are kept in one array per state;
 * sums over lineages are accumulated per lane, so results agree with Euler2ndOrder15 up to rounding.
 * Indicator pairs must come in the canonical order (a, b != a) as produced by the dynamics,
 * otherwise (and for very few lineages) Euler2ndOrderFused is used.
 * Requires the JVM option --add-modules jdk.incubator.vector.
 *
 * Generated code, do not edit by hand.
 */
public class Euler2ndOrderVec15 extends Euler2ndOrderFused {

	static final int N = 15;
	static final int NP = N * (N - 1);
	static final VectorSpecies<Double> SP = DoubleVector.SPECIES_PREFERRED;
	static final int W = SP.length();
	// migration routes (source, sink) in canonical order
	static final int[] SRC = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14};
	static final int[] SNK = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};

	final double[] rt = new double[NP], out = new double[N], cr = new double[N];
	final double[] sS = new double[N], sD = new double[N];
	double dLast, ddLast, maxDDD;
	double[] P0 = new double[0], P1 = new double[0], P2 = new double[0], P3 = new double[0], P4 = new double[0], P5 = new double[0], P6 = new double[0], P7 = new double[0], P8 = new double[0], P9 = new double[0], P10 = new double[0], P11 = new double[0], P12 = new double[0], P13 = new double[0], P14 = new double[0], D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, DD0, DD1, DD2, DD3, DD4, DD5, DD6, DD7, DD8, DD9, DD10, DD11, DD12, DD13, DD14;

	private void ensureCapacity(int size) {
		if (P0.length >= size) return;
		P0 = new double[size]; D0 = new double[size]; DD0 = new double[size];
		P1 = new double[size]; D1 = new double[size]; DD1 = new double[size];
		P2 = new double[size]; D2 = new double[size]; DD2 = new double[size];
		P3 = new double[size]; D3 = new double[size]; DD3 = new double[size];
		P4 = new double[size]; D4 = new double[size]; DD4 = new double[size];
		P5 = new double[size]; D5 = new double[size]; DD5 = new double[size];
		P6 = new double[size]; D6 = new double[size]; DD6 = new double[size];
		P7 = new double[size]; D7 = new double[size]; DD7 = new double[size];
		P8 = new double[size]; D8 = new double[size]; DD8 = new double[size];
		P9 = new double[size]; D9 = new double[size]; DD9 = new double[size];
		P10 = new double[size]; D10 = new double[size]; DD10 = new double[size];
		P11 = new double[size]; D11 = new double[size]; DD11 = new double[size];
		P12 = new double[size]; D12 = new double[size]; DD12 = new double[size];
		P13 = new double[size]; D13 = new double[size]; DD13 = new double[size];
		P14 = new double[size]; D14 = new double[size]; DD14 = new double[size];
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
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12, P13 = this.P13, P14 = this.P14;
		final double[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4, D5 = this.D5, D6 = this.D6, D7 = this.D7, D8 = this.D8, D9 = this.D9, D10 = this.D10, D11 = this.D11, D12 = this.D12, D13 = this.D13, D14 = this.D14;
		final double[] rt = this.rt, out = this.out;
		final double s0 = (double) sS[0], s1 = (double) sS[1], s2 = (double) sS[2], s3 = (double) sS[3], s4 = (double) sS[4], s5 = (double) sS[5], s6 = (double) sS[6], s7 = (double) sS[7], s8 = (double) sS[8], s9 = (double) sS[9], s10 = (double) sS[10], s11 = (double) sS[11], s12 = (double) sS[12], s13 = (double) sS[13], s14 = (double) sS[14];
		final double c0 = cr[0], c1 = cr[1], c2 = cr[2], c3 = cr[3], c4 = cr[4], c5 = cr[5], c6 = cr[6], c7 = cr[7], c8 = cr[8], c9 = cr[9], c10 = cr[10], c11 = cr[11], c12 = cr[12], c13 = cr[13], c14 = cr[14];
		DoubleVector accCoal = DoubleVector.zero(SP);
		DoubleVector aD0 = DoubleVector.zero(SP), aD1 = DoubleVector.zero(SP), aD2 = DoubleVector.zero(SP), aD3 = DoubleVector.zero(SP), aD4 = DoubleVector.zero(SP), aD5 = DoubleVector.zero(SP), aD6 = DoubleVector.zero(SP), aD7 = DoubleVector.zero(SP), aD8 = DoubleVector.zero(SP), aD9 = DoubleVector.zero(SP), aD10 = DoubleVector.zero(SP), aD11 = DoubleVector.zero(SP), aD12 = DoubleVector.zero(SP), aD13 = DoubleVector.zero(SP), aD14 = DoubleVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final DoubleVector p0 = DoubleVector.fromArray(SP, P0, i);
			final DoubleVector p1 = DoubleVector.fromArray(SP, P1, i);
			final DoubleVector p2 = DoubleVector.fromArray(SP, P2, i);
			final DoubleVector p3 = DoubleVector.fromArray(SP, P3, i);
			final DoubleVector p4 = DoubleVector.fromArray(SP, P4, i);
			final DoubleVector p5 = DoubleVector.fromArray(SP, P5, i);
			final DoubleVector p6 = DoubleVector.fromArray(SP, P6, i);
			final DoubleVector p7 = DoubleVector.fromArray(SP, P7, i);
			final DoubleVector p8 = DoubleVector.fromArray(SP, P8, i);
			final DoubleVector p9 = DoubleVector.fromArray(SP, P9, i);
			final DoubleVector p10 = DoubleVector.fromArray(SP, P10, i);
			final DoubleVector p11 = DoubleVector.fromArray(SP, P11, i);
			final DoubleVector p12 = DoubleVector.fromArray(SP, P12, i);
			final DoubleVector p13 = DoubleVector.fromArray(SP, P13, i);
			final DoubleVector p14 = DoubleVector.fromArray(SP, P14, i);
			final DoubleVector t0 = p0.neg().add(s0).mul(c0);
			final DoubleVector t1 = p1.neg().add(s1).mul(c1);
			final DoubleVector t2 = p2.neg().add(s2).mul(c2);
			final DoubleVector t3 = p3.neg().add(s3).mul(c3);
			final DoubleVector t4 = p4.neg().add(s4).mul(c4);
			final DoubleVector t5 = p5.neg().add(s5).mul(c5);
			final DoubleVector t6 = p6.neg().add(s6).mul(c6);
			final DoubleVector t7 = p7.neg().add(s7).mul(c7);
			final DoubleVector t8 = p8.neg().add(s8).mul(c8);
			final DoubleVector t9 = p9.neg().add(s9).mul(c9);
			final DoubleVector t10 = p10.neg().add(s10).mul(c10);
			final DoubleVector t11 = p11.neg().add(s11).mul(c11);
			final DoubleVector t12 = p12.neg().add(s12).mul(c12);
			final DoubleVector t13 = p13.neg().add(s13).mul(c13);
			final DoubleVector t14 = p14.neg().add(s14).mul(c14);
			final DoubleVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2)).add(p3.mul(t3)).add(p4.mul(t4)).add(p5.mul(t5)).add(p6.mul(t6)).add(p7.mul(t7)).add(p8.mul(t8)).add(p9.mul(t9)).add(p10.mul(t10)).add(p11.mul(t11)).add(p12.mul(t12)).add(p13.mul(t13)).add(p14.mul(t14));
			accCoal = accCoal.add(sc);
			final DoubleVector d0 = p0.mul(sc.sub(t0)).add(p1.mul(rt[14])).add(p2.mul(rt[28])).add(p3.mul(rt[42])).add(p4.mul(rt[56])).add(p5.mul(rt[70])).add(p6.mul(rt[84])).add(p7.mul(rt[98])).add(p8.mul(rt[112])).add(p9.mul(rt[126])).add(p10.mul(rt[140])).add(p11.mul(rt[154])).add(p12.mul(rt[168])).add(p13.mul(rt[182])).add(p14.mul(rt[196])).sub(p0.mul(out[0]));
			d0.intoArray(D0, i); aD0 = aD0.add(d0);
			final DoubleVector d1 = p1.mul(sc.sub(t1)).add(p0.mul(rt[0])).add(p2.mul(rt[29])).add(p3.mul(rt[43])).add(p4.mul(rt[57])).add(p5.mul(rt[71])).add(p6.mul(rt[85])).add(p7.mul(rt[99])).add(p8.mul(rt[113])).add(p9.mul(rt[127])).add(p10.mul(rt[141])).add(p11.mul(rt[155])).add(p12.mul(rt[169])).add(p13.mul(rt[183])).add(p14.mul(rt[197])).sub(p1.mul(out[1]));
			d1.intoArray(D1, i); aD1 = aD1.add(d1);
			final DoubleVector d2 = p2.mul(sc.sub(t2)).add(p0.mul(rt[1])).add(p1.mul(rt[15])).add(p3.mul(rt[44])).add(p4.mul(rt[58])).add(p5.mul(rt[72])).add(p6.mul(rt[86])).add(p7.mul(rt[100])).add(p8.mul(rt[114])).add(p9.mul(rt[128])).add(p10.mul(rt[142])).add(p11.mul(rt[156])).add(p12.mul(rt[170])).add(p13.mul(rt[184])).add(p14.mul(rt[198])).sub(p2.mul(out[2]));
			d2.intoArray(D2, i); aD2 = aD2.add(d2);
			final DoubleVector d3 = p3.mul(sc.sub(t3)).add(p0.mul(rt[2])).add(p1.mul(rt[16])).add(p2.mul(rt[30])).add(p4.mul(rt[59])).add(p5.mul(rt[73])).add(p6.mul(rt[87])).add(p7.mul(rt[101])).add(p8.mul(rt[115])).add(p9.mul(rt[129])).add(p10.mul(rt[143])).add(p11.mul(rt[157])).add(p12.mul(rt[171])).add(p13.mul(rt[185])).add(p14.mul(rt[199])).sub(p3.mul(out[3]));
			d3.intoArray(D3, i); aD3 = aD3.add(d3);
			final DoubleVector d4 = p4.mul(sc.sub(t4)).add(p0.mul(rt[3])).add(p1.mul(rt[17])).add(p2.mul(rt[31])).add(p3.mul(rt[45])).add(p5.mul(rt[74])).add(p6.mul(rt[88])).add(p7.mul(rt[102])).add(p8.mul(rt[116])).add(p9.mul(rt[130])).add(p10.mul(rt[144])).add(p11.mul(rt[158])).add(p12.mul(rt[172])).add(p13.mul(rt[186])).add(p14.mul(rt[200])).sub(p4.mul(out[4]));
			d4.intoArray(D4, i); aD4 = aD4.add(d4);
			final DoubleVector d5 = p5.mul(sc.sub(t5)).add(p0.mul(rt[4])).add(p1.mul(rt[18])).add(p2.mul(rt[32])).add(p3.mul(rt[46])).add(p4.mul(rt[60])).add(p6.mul(rt[89])).add(p7.mul(rt[103])).add(p8.mul(rt[117])).add(p9.mul(rt[131])).add(p10.mul(rt[145])).add(p11.mul(rt[159])).add(p12.mul(rt[173])).add(p13.mul(rt[187])).add(p14.mul(rt[201])).sub(p5.mul(out[5]));
			d5.intoArray(D5, i); aD5 = aD5.add(d5);
			final DoubleVector d6 = p6.mul(sc.sub(t6)).add(p0.mul(rt[5])).add(p1.mul(rt[19])).add(p2.mul(rt[33])).add(p3.mul(rt[47])).add(p4.mul(rt[61])).add(p5.mul(rt[75])).add(p7.mul(rt[104])).add(p8.mul(rt[118])).add(p9.mul(rt[132])).add(p10.mul(rt[146])).add(p11.mul(rt[160])).add(p12.mul(rt[174])).add(p13.mul(rt[188])).add(p14.mul(rt[202])).sub(p6.mul(out[6]));
			d6.intoArray(D6, i); aD6 = aD6.add(d6);
			final DoubleVector d7 = p7.mul(sc.sub(t7)).add(p0.mul(rt[6])).add(p1.mul(rt[20])).add(p2.mul(rt[34])).add(p3.mul(rt[48])).add(p4.mul(rt[62])).add(p5.mul(rt[76])).add(p6.mul(rt[90])).add(p8.mul(rt[119])).add(p9.mul(rt[133])).add(p10.mul(rt[147])).add(p11.mul(rt[161])).add(p12.mul(rt[175])).add(p13.mul(rt[189])).add(p14.mul(rt[203])).sub(p7.mul(out[7]));
			d7.intoArray(D7, i); aD7 = aD7.add(d7);
			final DoubleVector d8 = p8.mul(sc.sub(t8)).add(p0.mul(rt[7])).add(p1.mul(rt[21])).add(p2.mul(rt[35])).add(p3.mul(rt[49])).add(p4.mul(rt[63])).add(p5.mul(rt[77])).add(p6.mul(rt[91])).add(p7.mul(rt[105])).add(p9.mul(rt[134])).add(p10.mul(rt[148])).add(p11.mul(rt[162])).add(p12.mul(rt[176])).add(p13.mul(rt[190])).add(p14.mul(rt[204])).sub(p8.mul(out[8]));
			d8.intoArray(D8, i); aD8 = aD8.add(d8);
			final DoubleVector d9 = p9.mul(sc.sub(t9)).add(p0.mul(rt[8])).add(p1.mul(rt[22])).add(p2.mul(rt[36])).add(p3.mul(rt[50])).add(p4.mul(rt[64])).add(p5.mul(rt[78])).add(p6.mul(rt[92])).add(p7.mul(rt[106])).add(p8.mul(rt[120])).add(p10.mul(rt[149])).add(p11.mul(rt[163])).add(p12.mul(rt[177])).add(p13.mul(rt[191])).add(p14.mul(rt[205])).sub(p9.mul(out[9]));
			d9.intoArray(D9, i); aD9 = aD9.add(d9);
			final DoubleVector d10 = p10.mul(sc.sub(t10)).add(p0.mul(rt[9])).add(p1.mul(rt[23])).add(p2.mul(rt[37])).add(p3.mul(rt[51])).add(p4.mul(rt[65])).add(p5.mul(rt[79])).add(p6.mul(rt[93])).add(p7.mul(rt[107])).add(p8.mul(rt[121])).add(p9.mul(rt[135])).add(p11.mul(rt[164])).add(p12.mul(rt[178])).add(p13.mul(rt[192])).add(p14.mul(rt[206])).sub(p10.mul(out[10]));
			d10.intoArray(D10, i); aD10 = aD10.add(d10);
			final DoubleVector d11 = p11.mul(sc.sub(t11)).add(p0.mul(rt[10])).add(p1.mul(rt[24])).add(p2.mul(rt[38])).add(p3.mul(rt[52])).add(p4.mul(rt[66])).add(p5.mul(rt[80])).add(p6.mul(rt[94])).add(p7.mul(rt[108])).add(p8.mul(rt[122])).add(p9.mul(rt[136])).add(p10.mul(rt[150])).add(p12.mul(rt[179])).add(p13.mul(rt[193])).add(p14.mul(rt[207])).sub(p11.mul(out[11]));
			d11.intoArray(D11, i); aD11 = aD11.add(d11);
			final DoubleVector d12 = p12.mul(sc.sub(t12)).add(p0.mul(rt[11])).add(p1.mul(rt[25])).add(p2.mul(rt[39])).add(p3.mul(rt[53])).add(p4.mul(rt[67])).add(p5.mul(rt[81])).add(p6.mul(rt[95])).add(p7.mul(rt[109])).add(p8.mul(rt[123])).add(p9.mul(rt[137])).add(p10.mul(rt[151])).add(p11.mul(rt[165])).add(p13.mul(rt[194])).add(p14.mul(rt[208])).sub(p12.mul(out[12]));
			d12.intoArray(D12, i); aD12 = aD12.add(d12);
			final DoubleVector d13 = p13.mul(sc.sub(t13)).add(p0.mul(rt[12])).add(p1.mul(rt[26])).add(p2.mul(rt[40])).add(p3.mul(rt[54])).add(p4.mul(rt[68])).add(p5.mul(rt[82])).add(p6.mul(rt[96])).add(p7.mul(rt[110])).add(p8.mul(rt[124])).add(p9.mul(rt[138])).add(p10.mul(rt[152])).add(p11.mul(rt[166])).add(p12.mul(rt[180])).add(p14.mul(rt[209])).sub(p13.mul(out[13]));
			d13.intoArray(D13, i); aD13 = aD13.add(d13);
			final DoubleVector d14 = p14.mul(sc.sub(t14)).add(p0.mul(rt[13])).add(p1.mul(rt[27])).add(p2.mul(rt[41])).add(p3.mul(rt[55])).add(p4.mul(rt[69])).add(p5.mul(rt[83])).add(p6.mul(rt[97])).add(p7.mul(rt[111])).add(p8.mul(rt[125])).add(p9.mul(rt[139])).add(p10.mul(rt[153])).add(p11.mul(rt[167])).add(p12.mul(rt[181])).add(p13.mul(rt[195])).sub(p14.mul(out[14]));
			d14.intoArray(D14, i); aD14 = aD14.add(d14);
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
		sD[13] = aD13.reduceLanes(VectorOperators.ADD);
		sD[14] = aD14.reduceLanes(VectorOperators.ADD);
	}

	/** second derivative and the max of the approximate third derivative */
	private void pass2(final int Lpad) {
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12, P13 = this.P13, P14 = this.P14;
		final double[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4, D5 = this.D5, D6 = this.D6, D7 = this.D7, D8 = this.D8, D9 = this.D9, D10 = this.D10, D11 = this.D11, D12 = this.D12, D13 = this.D13, D14 = this.D14;
		final double[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4, DD5 = this.DD5, DD6 = this.DD6, DD7 = this.DD7, DD8 = this.DD8, DD9 = this.DD9, DD10 = this.DD10, DD11 = this.DD11, DD12 = this.DD12, DD13 = this.DD13, DD14 = this.DD14;
		final double[] rt = this.rt, out = this.out;
		final double s0 = (double) sS[0], s1 = (double) sS[1], s2 = (double) sS[2], s3 = (double) sS[3], s4 = (double) sS[4], s5 = (double) sS[5], s6 = (double) sS[6], s7 = (double) sS[7], s8 = (double) sS[8], s9 = (double) sS[9], s10 = (double) sS[10], s11 = (double) sS[11], s12 = (double) sS[12], s13 = (double) sS[13], s14 = (double) sS[14];
		final double c0 = cr[0], c1 = cr[1], c2 = cr[2], c3 = cr[3], c4 = cr[4], c5 = cr[5], c6 = cr[6], c7 = cr[7], c8 = cr[8], c9 = cr[9], c10 = cr[10], c11 = cr[11], c12 = cr[12], c13 = cr[13], c14 = cr[14];
		final double q0 = (double) sD[0], q1 = (double) sD[1], q2 = (double) sD[2], q3 = (double) sD[3], q4 = (double) sD[4], q5 = (double) sD[5], q6 = (double) sD[6], q7 = (double) sD[7], q8 = (double) sD[8], q9 = (double) sD[9], q10 = (double) sD[10], q11 = (double) sD[11], q12 = (double) sD[12], q13 = (double) sD[13], q14 = (double) sD[14];
		DoubleVector accPc = DoubleVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final DoubleVector p0 = DoubleVector.fromArray(SP, P0, i);
			final DoubleVector p1 = DoubleVector.fromArray(SP, P1, i);
			final DoubleVector p2 = DoubleVector.fromArray(SP, P2, i);
			final DoubleVector p3 = DoubleVector.fromArray(SP, P3, i);
			final DoubleVector p4 = DoubleVector.fromArray(SP, P4, i);
			final DoubleVector p5 = DoubleVector.fromArray(SP, P5, i);
			final DoubleVector p6 = DoubleVector.fromArray(SP, P6, i);
			final DoubleVector p7 = DoubleVector.fromArray(SP, P7, i);
			final DoubleVector p8 = DoubleVector.fromArray(SP, P8, i);
			final DoubleVector p9 = DoubleVector.fromArray(SP, P9, i);
			final DoubleVector p10 = DoubleVector.fromArray(SP, P10, i);
			final DoubleVector p11 = DoubleVector.fromArray(SP, P11, i);
			final DoubleVector p12 = DoubleVector.fromArray(SP, P12, i);
			final DoubleVector p13 = DoubleVector.fromArray(SP, P13, i);
			final DoubleVector p14 = DoubleVector.fromArray(SP, P14, i);
			final DoubleVector d0 = DoubleVector.fromArray(SP, D0, i);
			final DoubleVector d1 = DoubleVector.fromArray(SP, D1, i);
			final DoubleVector d2 = DoubleVector.fromArray(SP, D2, i);
			final DoubleVector d3 = DoubleVector.fromArray(SP, D3, i);
			final DoubleVector d4 = DoubleVector.fromArray(SP, D4, i);
			final DoubleVector d5 = DoubleVector.fromArray(SP, D5, i);
			final DoubleVector d6 = DoubleVector.fromArray(SP, D6, i);
			final DoubleVector d7 = DoubleVector.fromArray(SP, D7, i);
			final DoubleVector d8 = DoubleVector.fromArray(SP, D8, i);
			final DoubleVector d9 = DoubleVector.fromArray(SP, D9, i);
			final DoubleVector d10 = DoubleVector.fromArray(SP, D10, i);
			final DoubleVector d11 = DoubleVector.fromArray(SP, D11, i);
			final DoubleVector d12 = DoubleVector.fromArray(SP, D12, i);
			final DoubleVector d13 = DoubleVector.fromArray(SP, D13, i);
			final DoubleVector d14 = DoubleVector.fromArray(SP, D14, i);
			final DoubleVector t0 = p0.neg().add(s0).mul(c0);
			final DoubleVector t1 = p1.neg().add(s1).mul(c1);
			final DoubleVector t2 = p2.neg().add(s2).mul(c2);
			final DoubleVector t3 = p3.neg().add(s3).mul(c3);
			final DoubleVector t4 = p4.neg().add(s4).mul(c4);
			final DoubleVector t5 = p5.neg().add(s5).mul(c5);
			final DoubleVector t6 = p6.neg().add(s6).mul(c6);
			final DoubleVector t7 = p7.neg().add(s7).mul(c7);
			final DoubleVector t8 = p8.neg().add(s8).mul(c8);
			final DoubleVector t9 = p9.neg().add(s9).mul(c9);
			final DoubleVector t10 = p10.neg().add(s10).mul(c10);
			final DoubleVector t11 = p11.neg().add(s11).mul(c11);
			final DoubleVector t12 = p12.neg().add(s12).mul(c12);
			final DoubleVector t13 = p13.neg().add(s13).mul(c13);
			final DoubleVector t14 = p14.neg().add(s14).mul(c14);
			final DoubleVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2)).add(p3.mul(t3)).add(p4.mul(t4)).add(p5.mul(t5)).add(p6.mul(t6)).add(p7.mul(t7)).add(p8.mul(t8)).add(p9.mul(t9)).add(p10.mul(t10)).add(p11.mul(t11)).add(p12.mul(t12)).add(p13.mul(t13)).add(p14.mul(t14));
			final DoubleVector pc = d0.mul(p0.mul(-2.0).add(s0)).add(p0.mul(q0)).mul(c0).add(d1.mul(p1.mul(-2.0).add(s1)).add(p1.mul(q1)).mul(c1)).add(d2.mul(p2.mul(-2.0).add(s2)).add(p2.mul(q2)).mul(c2)).add(d3.mul(p3.mul(-2.0).add(s3)).add(p3.mul(q3)).mul(c3)).add(d4.mul(p4.mul(-2.0).add(s4)).add(p4.mul(q4)).mul(c4)).add(d5.mul(p5.mul(-2.0).add(s5)).add(p5.mul(q5)).mul(c5)).add(d6.mul(p6.mul(-2.0).add(s6)).add(p6.mul(q6)).mul(c6)).add(d7.mul(p7.mul(-2.0).add(s7)).add(p7.mul(q7)).mul(c7)).add(d8.mul(p8.mul(-2.0).add(s8)).add(p8.mul(q8)).mul(c8)).add(d9.mul(p9.mul(-2.0).add(s9)).add(p9.mul(q9)).mul(c9)).add(d10.mul(p10.mul(-2.0).add(s10)).add(p10.mul(q10)).mul(c10)).add(d11.mul(p11.mul(-2.0).add(s11)).add(p11.mul(q11)).mul(c11)).add(d12.mul(p12.mul(-2.0).add(s12)).add(p12.mul(q12)).mul(c12)).add(d13.mul(p13.mul(-2.0).add(s13)).add(p13.mul(q13)).mul(c13)).add(d14.mul(p14.mul(-2.0).add(s14)).add(p14.mul(q14)).mul(c14));
			accPc = accPc.add(pc);
			final DoubleVector dd0 = sc.sub(t0).mul(d0).add(p0.mul(pc.sub(d0.neg().add(q0).mul(c0)))).add(d1.mul(rt[14])).add(d2.mul(rt[28])).add(d3.mul(rt[42])).add(d4.mul(rt[56])).add(d5.mul(rt[70])).add(d6.mul(rt[84])).add(d7.mul(rt[98])).add(d8.mul(rt[112])).add(d9.mul(rt[126])).add(d10.mul(rt[140])).add(d11.mul(rt[154])).add(d12.mul(rt[168])).add(d13.mul(rt[182])).add(d14.mul(rt[196])).sub(d0.mul(out[0]));
			dd0.intoArray(DD0, i);
			final DoubleVector dd1 = sc.sub(t1).mul(d1).add(p1.mul(pc.sub(d1.neg().add(q1).mul(c1)))).add(d0.mul(rt[0])).add(d2.mul(rt[29])).add(d3.mul(rt[43])).add(d4.mul(rt[57])).add(d5.mul(rt[71])).add(d6.mul(rt[85])).add(d7.mul(rt[99])).add(d8.mul(rt[113])).add(d9.mul(rt[127])).add(d10.mul(rt[141])).add(d11.mul(rt[155])).add(d12.mul(rt[169])).add(d13.mul(rt[183])).add(d14.mul(rt[197])).sub(d1.mul(out[1]));
			dd1.intoArray(DD1, i);
			final DoubleVector dd2 = sc.sub(t2).mul(d2).add(p2.mul(pc.sub(d2.neg().add(q2).mul(c2)))).add(d0.mul(rt[1])).add(d1.mul(rt[15])).add(d3.mul(rt[44])).add(d4.mul(rt[58])).add(d5.mul(rt[72])).add(d6.mul(rt[86])).add(d7.mul(rt[100])).add(d8.mul(rt[114])).add(d9.mul(rt[128])).add(d10.mul(rt[142])).add(d11.mul(rt[156])).add(d12.mul(rt[170])).add(d13.mul(rt[184])).add(d14.mul(rt[198])).sub(d2.mul(out[2]));
			dd2.intoArray(DD2, i);
			final DoubleVector dd3 = sc.sub(t3).mul(d3).add(p3.mul(pc.sub(d3.neg().add(q3).mul(c3)))).add(d0.mul(rt[2])).add(d1.mul(rt[16])).add(d2.mul(rt[30])).add(d4.mul(rt[59])).add(d5.mul(rt[73])).add(d6.mul(rt[87])).add(d7.mul(rt[101])).add(d8.mul(rt[115])).add(d9.mul(rt[129])).add(d10.mul(rt[143])).add(d11.mul(rt[157])).add(d12.mul(rt[171])).add(d13.mul(rt[185])).add(d14.mul(rt[199])).sub(d3.mul(out[3]));
			dd3.intoArray(DD3, i);
			final DoubleVector dd4 = sc.sub(t4).mul(d4).add(p4.mul(pc.sub(d4.neg().add(q4).mul(c4)))).add(d0.mul(rt[3])).add(d1.mul(rt[17])).add(d2.mul(rt[31])).add(d3.mul(rt[45])).add(d5.mul(rt[74])).add(d6.mul(rt[88])).add(d7.mul(rt[102])).add(d8.mul(rt[116])).add(d9.mul(rt[130])).add(d10.mul(rt[144])).add(d11.mul(rt[158])).add(d12.mul(rt[172])).add(d13.mul(rt[186])).add(d14.mul(rt[200])).sub(d4.mul(out[4]));
			dd4.intoArray(DD4, i);
			final DoubleVector dd5 = sc.sub(t5).mul(d5).add(p5.mul(pc.sub(d5.neg().add(q5).mul(c5)))).add(d0.mul(rt[4])).add(d1.mul(rt[18])).add(d2.mul(rt[32])).add(d3.mul(rt[46])).add(d4.mul(rt[60])).add(d6.mul(rt[89])).add(d7.mul(rt[103])).add(d8.mul(rt[117])).add(d9.mul(rt[131])).add(d10.mul(rt[145])).add(d11.mul(rt[159])).add(d12.mul(rt[173])).add(d13.mul(rt[187])).add(d14.mul(rt[201])).sub(d5.mul(out[5]));
			dd5.intoArray(DD5, i);
			final DoubleVector dd6 = sc.sub(t6).mul(d6).add(p6.mul(pc.sub(d6.neg().add(q6).mul(c6)))).add(d0.mul(rt[5])).add(d1.mul(rt[19])).add(d2.mul(rt[33])).add(d3.mul(rt[47])).add(d4.mul(rt[61])).add(d5.mul(rt[75])).add(d7.mul(rt[104])).add(d8.mul(rt[118])).add(d9.mul(rt[132])).add(d10.mul(rt[146])).add(d11.mul(rt[160])).add(d12.mul(rt[174])).add(d13.mul(rt[188])).add(d14.mul(rt[202])).sub(d6.mul(out[6]));
			dd6.intoArray(DD6, i);
			final DoubleVector dd7 = sc.sub(t7).mul(d7).add(p7.mul(pc.sub(d7.neg().add(q7).mul(c7)))).add(d0.mul(rt[6])).add(d1.mul(rt[20])).add(d2.mul(rt[34])).add(d3.mul(rt[48])).add(d4.mul(rt[62])).add(d5.mul(rt[76])).add(d6.mul(rt[90])).add(d8.mul(rt[119])).add(d9.mul(rt[133])).add(d10.mul(rt[147])).add(d11.mul(rt[161])).add(d12.mul(rt[175])).add(d13.mul(rt[189])).add(d14.mul(rt[203])).sub(d7.mul(out[7]));
			dd7.intoArray(DD7, i);
			final DoubleVector dd8 = sc.sub(t8).mul(d8).add(p8.mul(pc.sub(d8.neg().add(q8).mul(c8)))).add(d0.mul(rt[7])).add(d1.mul(rt[21])).add(d2.mul(rt[35])).add(d3.mul(rt[49])).add(d4.mul(rt[63])).add(d5.mul(rt[77])).add(d6.mul(rt[91])).add(d7.mul(rt[105])).add(d9.mul(rt[134])).add(d10.mul(rt[148])).add(d11.mul(rt[162])).add(d12.mul(rt[176])).add(d13.mul(rt[190])).add(d14.mul(rt[204])).sub(d8.mul(out[8]));
			dd8.intoArray(DD8, i);
			final DoubleVector dd9 = sc.sub(t9).mul(d9).add(p9.mul(pc.sub(d9.neg().add(q9).mul(c9)))).add(d0.mul(rt[8])).add(d1.mul(rt[22])).add(d2.mul(rt[36])).add(d3.mul(rt[50])).add(d4.mul(rt[64])).add(d5.mul(rt[78])).add(d6.mul(rt[92])).add(d7.mul(rt[106])).add(d8.mul(rt[120])).add(d10.mul(rt[149])).add(d11.mul(rt[163])).add(d12.mul(rt[177])).add(d13.mul(rt[191])).add(d14.mul(rt[205])).sub(d9.mul(out[9]));
			dd9.intoArray(DD9, i);
			final DoubleVector dd10 = sc.sub(t10).mul(d10).add(p10.mul(pc.sub(d10.neg().add(q10).mul(c10)))).add(d0.mul(rt[9])).add(d1.mul(rt[23])).add(d2.mul(rt[37])).add(d3.mul(rt[51])).add(d4.mul(rt[65])).add(d5.mul(rt[79])).add(d6.mul(rt[93])).add(d7.mul(rt[107])).add(d8.mul(rt[121])).add(d9.mul(rt[135])).add(d11.mul(rt[164])).add(d12.mul(rt[178])).add(d13.mul(rt[192])).add(d14.mul(rt[206])).sub(d10.mul(out[10]));
			dd10.intoArray(DD10, i);
			final DoubleVector dd11 = sc.sub(t11).mul(d11).add(p11.mul(pc.sub(d11.neg().add(q11).mul(c11)))).add(d0.mul(rt[10])).add(d1.mul(rt[24])).add(d2.mul(rt[38])).add(d3.mul(rt[52])).add(d4.mul(rt[66])).add(d5.mul(rt[80])).add(d6.mul(rt[94])).add(d7.mul(rt[108])).add(d8.mul(rt[122])).add(d9.mul(rt[136])).add(d10.mul(rt[150])).add(d12.mul(rt[179])).add(d13.mul(rt[193])).add(d14.mul(rt[207])).sub(d11.mul(out[11]));
			dd11.intoArray(DD11, i);
			final DoubleVector dd12 = sc.sub(t12).mul(d12).add(p12.mul(pc.sub(d12.neg().add(q12).mul(c12)))).add(d0.mul(rt[11])).add(d1.mul(rt[25])).add(d2.mul(rt[39])).add(d3.mul(rt[53])).add(d4.mul(rt[67])).add(d5.mul(rt[81])).add(d6.mul(rt[95])).add(d7.mul(rt[109])).add(d8.mul(rt[123])).add(d9.mul(rt[137])).add(d10.mul(rt[151])).add(d11.mul(rt[165])).add(d13.mul(rt[194])).add(d14.mul(rt[208])).sub(d12.mul(out[12]));
			dd12.intoArray(DD12, i);
			final DoubleVector dd13 = sc.sub(t13).mul(d13).add(p13.mul(pc.sub(d13.neg().add(q13).mul(c13)))).add(d0.mul(rt[12])).add(d1.mul(rt[26])).add(d2.mul(rt[40])).add(d3.mul(rt[54])).add(d4.mul(rt[68])).add(d5.mul(rt[82])).add(d6.mul(rt[96])).add(d7.mul(rt[110])).add(d8.mul(rt[124])).add(d9.mul(rt[138])).add(d10.mul(rt[152])).add(d11.mul(rt[166])).add(d12.mul(rt[180])).add(d14.mul(rt[209])).sub(d13.mul(out[13]));
			dd13.intoArray(DD13, i);
			final DoubleVector dd14 = sc.sub(t14).mul(d14).add(p14.mul(pc.sub(d14.neg().add(q14).mul(c14)))).add(d0.mul(rt[13])).add(d1.mul(rt[27])).add(d2.mul(rt[41])).add(d3.mul(rt[55])).add(d4.mul(rt[69])).add(d5.mul(rt[83])).add(d6.mul(rt[97])).add(d7.mul(rt[111])).add(d8.mul(rt[125])).add(d9.mul(rt[139])).add(d10.mul(rt[153])).add(d11.mul(rt[167])).add(d12.mul(rt[181])).add(d13.mul(rt[195])).sub(d14.mul(out[14]));
			dd14.intoArray(DD14, i);
		}
		ddLast = -accPc.reduceLanes(VectorOperators.ADD) / 2;
	}

	/** max of the approximate third derivative */
	private void pass2b(final int Lpad) {
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12, P13 = this.P13, P14 = this.P14;
		final double[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4, DD5 = this.DD5, DD6 = this.DD6, DD7 = this.DD7, DD8 = this.DD8, DD9 = this.DD9, DD10 = this.DD10, DD11 = this.DD11, DD12 = this.DD12, DD13 = this.DD13, DD14 = this.DD14;
		final double[] rt = this.rt, out = this.out;
		final double s0 = (double) sS[0], s1 = (double) sS[1], s2 = (double) sS[2], s3 = (double) sS[3], s4 = (double) sS[4], s5 = (double) sS[5], s6 = (double) sS[6], s7 = (double) sS[7], s8 = (double) sS[8], s9 = (double) sS[9], s10 = (double) sS[10], s11 = (double) sS[11], s12 = (double) sS[12], s13 = (double) sS[13], s14 = (double) sS[14];
		final double c0 = cr[0], c1 = cr[1], c2 = cr[2], c3 = cr[3], c4 = cr[4], c5 = cr[5], c6 = cr[6], c7 = cr[7], c8 = cr[8], c9 = cr[9], c10 = cr[10], c11 = cr[11], c12 = cr[12], c13 = cr[13], c14 = cr[14];
		DoubleVector accMax = DoubleVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final DoubleVector p0 = DoubleVector.fromArray(SP, P0, i);
			final DoubleVector p1 = DoubleVector.fromArray(SP, P1, i);
			final DoubleVector p2 = DoubleVector.fromArray(SP, P2, i);
			final DoubleVector p3 = DoubleVector.fromArray(SP, P3, i);
			final DoubleVector p4 = DoubleVector.fromArray(SP, P4, i);
			final DoubleVector p5 = DoubleVector.fromArray(SP, P5, i);
			final DoubleVector p6 = DoubleVector.fromArray(SP, P6, i);
			final DoubleVector p7 = DoubleVector.fromArray(SP, P7, i);
			final DoubleVector p8 = DoubleVector.fromArray(SP, P8, i);
			final DoubleVector p9 = DoubleVector.fromArray(SP, P9, i);
			final DoubleVector p10 = DoubleVector.fromArray(SP, P10, i);
			final DoubleVector p11 = DoubleVector.fromArray(SP, P11, i);
			final DoubleVector p12 = DoubleVector.fromArray(SP, P12, i);
			final DoubleVector p13 = DoubleVector.fromArray(SP, P13, i);
			final DoubleVector p14 = DoubleVector.fromArray(SP, P14, i);
			final DoubleVector dd0 = DoubleVector.fromArray(SP, DD0, i);
			final DoubleVector dd1 = DoubleVector.fromArray(SP, DD1, i);
			final DoubleVector dd2 = DoubleVector.fromArray(SP, DD2, i);
			final DoubleVector dd3 = DoubleVector.fromArray(SP, DD3, i);
			final DoubleVector dd4 = DoubleVector.fromArray(SP, DD4, i);
			final DoubleVector dd5 = DoubleVector.fromArray(SP, DD5, i);
			final DoubleVector dd6 = DoubleVector.fromArray(SP, DD6, i);
			final DoubleVector dd7 = DoubleVector.fromArray(SP, DD7, i);
			final DoubleVector dd8 = DoubleVector.fromArray(SP, DD8, i);
			final DoubleVector dd9 = DoubleVector.fromArray(SP, DD9, i);
			final DoubleVector dd10 = DoubleVector.fromArray(SP, DD10, i);
			final DoubleVector dd11 = DoubleVector.fromArray(SP, DD11, i);
			final DoubleVector dd12 = DoubleVector.fromArray(SP, DD12, i);
			final DoubleVector dd13 = DoubleVector.fromArray(SP, DD13, i);
			final DoubleVector dd14 = DoubleVector.fromArray(SP, DD14, i);
			final DoubleVector t0 = p0.neg().add(s0).mul(c0);
			final DoubleVector t1 = p1.neg().add(s1).mul(c1);
			final DoubleVector t2 = p2.neg().add(s2).mul(c2);
			final DoubleVector t3 = p3.neg().add(s3).mul(c3);
			final DoubleVector t4 = p4.neg().add(s4).mul(c4);
			final DoubleVector t5 = p5.neg().add(s5).mul(c5);
			final DoubleVector t6 = p6.neg().add(s6).mul(c6);
			final DoubleVector t7 = p7.neg().add(s7).mul(c7);
			final DoubleVector t8 = p8.neg().add(s8).mul(c8);
			final DoubleVector t9 = p9.neg().add(s9).mul(c9);
			final DoubleVector t10 = p10.neg().add(s10).mul(c10);
			final DoubleVector t11 = p11.neg().add(s11).mul(c11);
			final DoubleVector t12 = p12.neg().add(s12).mul(c12);
			final DoubleVector t13 = p13.neg().add(s13).mul(c13);
			final DoubleVector t14 = p14.neg().add(s14).mul(c14);
			final DoubleVector sc = p0.mul(t0).add(p1.mul(t1)).add(p2.mul(t2)).add(p3.mul(t3)).add(p4.mul(t4)).add(p5.mul(t5)).add(p6.mul(t6)).add(p7.mul(t7)).add(p8.mul(t8)).add(p9.mul(t9)).add(p10.mul(t10)).add(p11.mul(t11)).add(p12.mul(t12)).add(p13.mul(t13)).add(p14.mul(t14));
			accMax = accMax.max(sc.sub(t0).mul(dd0).add(dd1.mul(rt[14])).add(dd2.mul(rt[28])).add(dd3.mul(rt[42])).add(dd4.mul(rt[56])).add(dd5.mul(rt[70])).add(dd6.mul(rt[84])).add(dd7.mul(rt[98])).add(dd8.mul(rt[112])).add(dd9.mul(rt[126])).add(dd10.mul(rt[140])).add(dd11.mul(rt[154])).add(dd12.mul(rt[168])).add(dd13.mul(rt[182])).add(dd14.mul(rt[196])).sub(dd0.mul(out[0])).abs());
			accMax = accMax.max(sc.sub(t1).mul(dd1).add(dd0.mul(rt[0])).add(dd2.mul(rt[29])).add(dd3.mul(rt[43])).add(dd4.mul(rt[57])).add(dd5.mul(rt[71])).add(dd6.mul(rt[85])).add(dd7.mul(rt[99])).add(dd8.mul(rt[113])).add(dd9.mul(rt[127])).add(dd10.mul(rt[141])).add(dd11.mul(rt[155])).add(dd12.mul(rt[169])).add(dd13.mul(rt[183])).add(dd14.mul(rt[197])).sub(dd1.mul(out[1])).abs());
			accMax = accMax.max(sc.sub(t2).mul(dd2).add(dd0.mul(rt[1])).add(dd1.mul(rt[15])).add(dd3.mul(rt[44])).add(dd4.mul(rt[58])).add(dd5.mul(rt[72])).add(dd6.mul(rt[86])).add(dd7.mul(rt[100])).add(dd8.mul(rt[114])).add(dd9.mul(rt[128])).add(dd10.mul(rt[142])).add(dd11.mul(rt[156])).add(dd12.mul(rt[170])).add(dd13.mul(rt[184])).add(dd14.mul(rt[198])).sub(dd2.mul(out[2])).abs());
			accMax = accMax.max(sc.sub(t3).mul(dd3).add(dd0.mul(rt[2])).add(dd1.mul(rt[16])).add(dd2.mul(rt[30])).add(dd4.mul(rt[59])).add(dd5.mul(rt[73])).add(dd6.mul(rt[87])).add(dd7.mul(rt[101])).add(dd8.mul(rt[115])).add(dd9.mul(rt[129])).add(dd10.mul(rt[143])).add(dd11.mul(rt[157])).add(dd12.mul(rt[171])).add(dd13.mul(rt[185])).add(dd14.mul(rt[199])).sub(dd3.mul(out[3])).abs());
			accMax = accMax.max(sc.sub(t4).mul(dd4).add(dd0.mul(rt[3])).add(dd1.mul(rt[17])).add(dd2.mul(rt[31])).add(dd3.mul(rt[45])).add(dd5.mul(rt[74])).add(dd6.mul(rt[88])).add(dd7.mul(rt[102])).add(dd8.mul(rt[116])).add(dd9.mul(rt[130])).add(dd10.mul(rt[144])).add(dd11.mul(rt[158])).add(dd12.mul(rt[172])).add(dd13.mul(rt[186])).add(dd14.mul(rt[200])).sub(dd4.mul(out[4])).abs());
			accMax = accMax.max(sc.sub(t5).mul(dd5).add(dd0.mul(rt[4])).add(dd1.mul(rt[18])).add(dd2.mul(rt[32])).add(dd3.mul(rt[46])).add(dd4.mul(rt[60])).add(dd6.mul(rt[89])).add(dd7.mul(rt[103])).add(dd8.mul(rt[117])).add(dd9.mul(rt[131])).add(dd10.mul(rt[145])).add(dd11.mul(rt[159])).add(dd12.mul(rt[173])).add(dd13.mul(rt[187])).add(dd14.mul(rt[201])).sub(dd5.mul(out[5])).abs());
			accMax = accMax.max(sc.sub(t6).mul(dd6).add(dd0.mul(rt[5])).add(dd1.mul(rt[19])).add(dd2.mul(rt[33])).add(dd3.mul(rt[47])).add(dd4.mul(rt[61])).add(dd5.mul(rt[75])).add(dd7.mul(rt[104])).add(dd8.mul(rt[118])).add(dd9.mul(rt[132])).add(dd10.mul(rt[146])).add(dd11.mul(rt[160])).add(dd12.mul(rt[174])).add(dd13.mul(rt[188])).add(dd14.mul(rt[202])).sub(dd6.mul(out[6])).abs());
			accMax = accMax.max(sc.sub(t7).mul(dd7).add(dd0.mul(rt[6])).add(dd1.mul(rt[20])).add(dd2.mul(rt[34])).add(dd3.mul(rt[48])).add(dd4.mul(rt[62])).add(dd5.mul(rt[76])).add(dd6.mul(rt[90])).add(dd8.mul(rt[119])).add(dd9.mul(rt[133])).add(dd10.mul(rt[147])).add(dd11.mul(rt[161])).add(dd12.mul(rt[175])).add(dd13.mul(rt[189])).add(dd14.mul(rt[203])).sub(dd7.mul(out[7])).abs());
			accMax = accMax.max(sc.sub(t8).mul(dd8).add(dd0.mul(rt[7])).add(dd1.mul(rt[21])).add(dd2.mul(rt[35])).add(dd3.mul(rt[49])).add(dd4.mul(rt[63])).add(dd5.mul(rt[77])).add(dd6.mul(rt[91])).add(dd7.mul(rt[105])).add(dd9.mul(rt[134])).add(dd10.mul(rt[148])).add(dd11.mul(rt[162])).add(dd12.mul(rt[176])).add(dd13.mul(rt[190])).add(dd14.mul(rt[204])).sub(dd8.mul(out[8])).abs());
			accMax = accMax.max(sc.sub(t9).mul(dd9).add(dd0.mul(rt[8])).add(dd1.mul(rt[22])).add(dd2.mul(rt[36])).add(dd3.mul(rt[50])).add(dd4.mul(rt[64])).add(dd5.mul(rt[78])).add(dd6.mul(rt[92])).add(dd7.mul(rt[106])).add(dd8.mul(rt[120])).add(dd10.mul(rt[149])).add(dd11.mul(rt[163])).add(dd12.mul(rt[177])).add(dd13.mul(rt[191])).add(dd14.mul(rt[205])).sub(dd9.mul(out[9])).abs());
			accMax = accMax.max(sc.sub(t10).mul(dd10).add(dd0.mul(rt[9])).add(dd1.mul(rt[23])).add(dd2.mul(rt[37])).add(dd3.mul(rt[51])).add(dd4.mul(rt[65])).add(dd5.mul(rt[79])).add(dd6.mul(rt[93])).add(dd7.mul(rt[107])).add(dd8.mul(rt[121])).add(dd9.mul(rt[135])).add(dd11.mul(rt[164])).add(dd12.mul(rt[178])).add(dd13.mul(rt[192])).add(dd14.mul(rt[206])).sub(dd10.mul(out[10])).abs());
			accMax = accMax.max(sc.sub(t11).mul(dd11).add(dd0.mul(rt[10])).add(dd1.mul(rt[24])).add(dd2.mul(rt[38])).add(dd3.mul(rt[52])).add(dd4.mul(rt[66])).add(dd5.mul(rt[80])).add(dd6.mul(rt[94])).add(dd7.mul(rt[108])).add(dd8.mul(rt[122])).add(dd9.mul(rt[136])).add(dd10.mul(rt[150])).add(dd12.mul(rt[179])).add(dd13.mul(rt[193])).add(dd14.mul(rt[207])).sub(dd11.mul(out[11])).abs());
			accMax = accMax.max(sc.sub(t12).mul(dd12).add(dd0.mul(rt[11])).add(dd1.mul(rt[25])).add(dd2.mul(rt[39])).add(dd3.mul(rt[53])).add(dd4.mul(rt[67])).add(dd5.mul(rt[81])).add(dd6.mul(rt[95])).add(dd7.mul(rt[109])).add(dd8.mul(rt[123])).add(dd9.mul(rt[137])).add(dd10.mul(rt[151])).add(dd11.mul(rt[165])).add(dd13.mul(rt[194])).add(dd14.mul(rt[208])).sub(dd12.mul(out[12])).abs());
			accMax = accMax.max(sc.sub(t13).mul(dd13).add(dd0.mul(rt[12])).add(dd1.mul(rt[26])).add(dd2.mul(rt[40])).add(dd3.mul(rt[54])).add(dd4.mul(rt[68])).add(dd5.mul(rt[82])).add(dd6.mul(rt[96])).add(dd7.mul(rt[110])).add(dd8.mul(rt[124])).add(dd9.mul(rt[138])).add(dd10.mul(rt[152])).add(dd11.mul(rt[166])).add(dd12.mul(rt[180])).add(dd14.mul(rt[209])).sub(dd13.mul(out[13])).abs());
			accMax = accMax.max(sc.sub(t14).mul(dd14).add(dd0.mul(rt[13])).add(dd1.mul(rt[27])).add(dd2.mul(rt[41])).add(dd3.mul(rt[55])).add(dd4.mul(rt[69])).add(dd5.mul(rt[83])).add(dd6.mul(rt[97])).add(dd7.mul(rt[111])).add(dd8.mul(rt[125])).add(dd9.mul(rt[139])).add(dd10.mul(rt[153])).add(dd11.mul(rt[167])).add(dd12.mul(rt[181])).add(dd13.mul(rt[195])).sub(dd14.mul(out[14])).abs());
		}
		maxDDD = accMax.reduceLanes(VectorOperators.MAX);
	}

	/** true if any probability would leave [0,1] or change by more than 0.2 with this step */
	private boolean needsSearch(final int Lpad, final double ts, final double tss) {
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12, P13 = this.P13, P14 = this.P14;
		final double[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4, D5 = this.D5, D6 = this.D6, D7 = this.D7, D8 = this.D8, D9 = this.D9, D10 = this.D10, D11 = this.D11, D12 = this.D12, D13 = this.D13, D14 = this.D14;
		final double[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4, DD5 = this.DD5, DD6 = this.DD6, DD7 = this.DD7, DD8 = this.DD8, DD9 = this.DD9, DD10 = this.DD10, DD11 = this.DD11, DD12 = this.DD12, DD13 = this.DD13, DD14 = this.DD14;
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
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P5, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D5, i).mul(ts)).add(DoubleVector.fromArray(SP, DD5, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P6, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D6, i).mul(ts)).add(DoubleVector.fromArray(SP, DD6, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P7, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D7, i).mul(ts)).add(DoubleVector.fromArray(SP, DD7, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P8, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D8, i).mul(ts)).add(DoubleVector.fromArray(SP, DD8, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P9, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D9, i).mul(ts)).add(DoubleVector.fromArray(SP, DD9, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P10, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D10, i).mul(ts)).add(DoubleVector.fromArray(SP, DD10, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P11, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D11, i).mul(ts)).add(DoubleVector.fromArray(SP, DD11, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P12, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D12, i).mul(ts)).add(DoubleVector.fromArray(SP, DD12, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P13, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D13, i).mul(ts)).add(DoubleVector.fromArray(SP, DD13, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			{
				final DoubleVector p = DoubleVector.fromArray(SP, P14, i);
				final DoubleVector nv = p.add(DoubleVector.fromArray(SP, D14, i).mul(ts)).add(DoubleVector.fromArray(SP, DD14, i).mul(tss));
				m = m.or(nv.compare(VectorOperators.GT, 1.0)).or(nv.compare(VectorOperators.LT, 0.0)).or(nv.sub(p).abs().compare(VectorOperators.GT, 0.2));
			}
			if (m.anyTrue()) return true;
		}
		return false;
	}

	/** update and normalise; returns true if a probability became negative */
	private boolean pass4(final int Lpad, final int L, final double ts, final double tss) {
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12, P13 = this.P13, P14 = this.P14;
		final double[] D0 = this.D0, D1 = this.D1, D2 = this.D2, D3 = this.D3, D4 = this.D4, D5 = this.D5, D6 = this.D6, D7 = this.D7, D8 = this.D8, D9 = this.D9, D10 = this.D10, D11 = this.D11, D12 = this.D12, D13 = this.D13, D14 = this.D14;
		final double[] DD0 = this.DD0, DD1 = this.DD1, DD2 = this.DD2, DD3 = this.DD3, DD4 = this.DD4, DD5 = this.DD5, DD6 = this.DD6, DD7 = this.DD7, DD8 = this.DD8, DD9 = this.DD9, DD10 = this.DD10, DD11 = this.DD11, DD12 = this.DD12, DD13 = this.DD13, DD14 = this.DD14;
		final int lastBlock = Lpad - W;
		final VectorMask<Double> pad = SP.indexInRange(lastBlock, L).not();
		boolean failed = false;
		DoubleVector aS0 = DoubleVector.zero(SP), aS1 = DoubleVector.zero(SP), aS2 = DoubleVector.zero(SP), aS3 = DoubleVector.zero(SP), aS4 = DoubleVector.zero(SP), aS5 = DoubleVector.zero(SP), aS6 = DoubleVector.zero(SP), aS7 = DoubleVector.zero(SP), aS8 = DoubleVector.zero(SP), aS9 = DoubleVector.zero(SP), aS10 = DoubleVector.zero(SP), aS11 = DoubleVector.zero(SP), aS12 = DoubleVector.zero(SP), aS13 = DoubleVector.zero(SP), aS14 = DoubleVector.zero(SP);
		for (int i = 0; i < Lpad; i += W) {
			final DoubleVector x0 = DoubleVector.fromArray(SP, P0, i).add(DoubleVector.fromArray(SP, D0, i).mul(ts).add(DoubleVector.fromArray(SP, DD0, i).mul(tss)));
			final DoubleVector x1 = DoubleVector.fromArray(SP, P1, i).add(DoubleVector.fromArray(SP, D1, i).mul(ts).add(DoubleVector.fromArray(SP, DD1, i).mul(tss)));
			final DoubleVector x2 = DoubleVector.fromArray(SP, P2, i).add(DoubleVector.fromArray(SP, D2, i).mul(ts).add(DoubleVector.fromArray(SP, DD2, i).mul(tss)));
			final DoubleVector x3 = DoubleVector.fromArray(SP, P3, i).add(DoubleVector.fromArray(SP, D3, i).mul(ts).add(DoubleVector.fromArray(SP, DD3, i).mul(tss)));
			final DoubleVector x4 = DoubleVector.fromArray(SP, P4, i).add(DoubleVector.fromArray(SP, D4, i).mul(ts).add(DoubleVector.fromArray(SP, DD4, i).mul(tss)));
			final DoubleVector x5 = DoubleVector.fromArray(SP, P5, i).add(DoubleVector.fromArray(SP, D5, i).mul(ts).add(DoubleVector.fromArray(SP, DD5, i).mul(tss)));
			final DoubleVector x6 = DoubleVector.fromArray(SP, P6, i).add(DoubleVector.fromArray(SP, D6, i).mul(ts).add(DoubleVector.fromArray(SP, DD6, i).mul(tss)));
			final DoubleVector x7 = DoubleVector.fromArray(SP, P7, i).add(DoubleVector.fromArray(SP, D7, i).mul(ts).add(DoubleVector.fromArray(SP, DD7, i).mul(tss)));
			final DoubleVector x8 = DoubleVector.fromArray(SP, P8, i).add(DoubleVector.fromArray(SP, D8, i).mul(ts).add(DoubleVector.fromArray(SP, DD8, i).mul(tss)));
			final DoubleVector x9 = DoubleVector.fromArray(SP, P9, i).add(DoubleVector.fromArray(SP, D9, i).mul(ts).add(DoubleVector.fromArray(SP, DD9, i).mul(tss)));
			final DoubleVector x10 = DoubleVector.fromArray(SP, P10, i).add(DoubleVector.fromArray(SP, D10, i).mul(ts).add(DoubleVector.fromArray(SP, DD10, i).mul(tss)));
			final DoubleVector x11 = DoubleVector.fromArray(SP, P11, i).add(DoubleVector.fromArray(SP, D11, i).mul(ts).add(DoubleVector.fromArray(SP, DD11, i).mul(tss)));
			final DoubleVector x12 = DoubleVector.fromArray(SP, P12, i).add(DoubleVector.fromArray(SP, D12, i).mul(ts).add(DoubleVector.fromArray(SP, DD12, i).mul(tss)));
			final DoubleVector x13 = DoubleVector.fromArray(SP, P13, i).add(DoubleVector.fromArray(SP, D13, i).mul(ts).add(DoubleVector.fromArray(SP, DD13, i).mul(tss)));
			final DoubleVector x14 = DoubleVector.fromArray(SP, P14, i).add(DoubleVector.fromArray(SP, D14, i).mul(ts).add(DoubleVector.fromArray(SP, DD14, i).mul(tss)));
			final VectorMask<Double> negative = x0.compare(VectorOperators.LT, 0.0).or(x1.compare(VectorOperators.LT, 0.0)).or(x2.compare(VectorOperators.LT, 0.0)).or(x3.compare(VectorOperators.LT, 0.0)).or(x4.compare(VectorOperators.LT, 0.0)).or(x5.compare(VectorOperators.LT, 0.0)).or(x6.compare(VectorOperators.LT, 0.0)).or(x7.compare(VectorOperators.LT, 0.0)).or(x8.compare(VectorOperators.LT, 0.0)).or(x9.compare(VectorOperators.LT, 0.0)).or(x10.compare(VectorOperators.LT, 0.0)).or(x11.compare(VectorOperators.LT, 0.0)).or(x12.compare(VectorOperators.LT, 0.0)).or(x13.compare(VectorOperators.LT, 0.0)).or(x14.compare(VectorOperators.LT, 0.0));
			final DoubleVector sum = x0.add(x1).add(x2).add(x3).add(x4).add(x5).add(x6).add(x7).add(x8).add(x9).add(x10).add(x11).add(x12).add(x13).add(x14);
			final boolean anyNegative = negative.anyTrue();
			if (anyNegative) failed = true;
			DoubleVector y0 = x0.div(sum);
			DoubleVector y1 = x1.div(sum);
			DoubleVector y2 = x2.div(sum);
			DoubleVector y3 = x3.div(sum);
			DoubleVector y4 = x4.div(sum);
			DoubleVector y5 = x5.div(sum);
			DoubleVector y6 = x6.div(sum);
			DoubleVector y7 = x7.div(sum);
			DoubleVector y8 = x8.div(sum);
			DoubleVector y9 = x9.div(sum);
			DoubleVector y10 = x10.div(sum);
			DoubleVector y11 = x11.div(sum);
			DoubleVector y12 = x12.div(sum);
			DoubleVector y13 = x13.div(sum);
			DoubleVector y14 = x14.div(sum);
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
				y13 = y13.blend(x13, negative);
				y14 = y14.blend(x14, negative);
			}
			if (i == lastBlock) {
				// padding lineages stay at zero (0/0 would give NaN)
				y0 = y0.blend(0.0, pad);
				y1 = y1.blend(0.0, pad);
				y2 = y2.blend(0.0, pad);
				y3 = y3.blend(0.0, pad);
				y4 = y4.blend(0.0, pad);
				y5 = y5.blend(0.0, pad);
				y6 = y6.blend(0.0, pad);
				y7 = y7.blend(0.0, pad);
				y8 = y8.blend(0.0, pad);
				y9 = y9.blend(0.0, pad);
				y10 = y10.blend(0.0, pad);
				y11 = y11.blend(0.0, pad);
				y12 = y12.blend(0.0, pad);
				y13 = y13.blend(0.0, pad);
				y14 = y14.blend(0.0, pad);
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
			y13.intoArray(P13, i); aS13 = aS13.add(y13);
			y14.intoArray(P14, i); aS14 = aS14.add(y14);
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
		sS[13] = aS13.reduceLanes(VectorOperators.ADD);
		sS[14] = aS14.reduceLanes(VectorOperators.ADD);
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
		final double[] P0 = this.P0, P1 = this.P1, P2 = this.P2, P3 = this.P3, P4 = this.P4, P5 = this.P5, P6 = this.P6, P7 = this.P7, P8 = this.P8, P9 = this.P9, P10 = this.P10, P11 = this.P11, P12 = this.P12, P13 = this.P13, P14 = this.P14;
		for (int i = 0, k = 0; i < L; i++, k += N) {
			P0[i] = (double) p[k + 0];
			P1[i] = (double) p[k + 1];
			P2[i] = (double) p[k + 2];
			P3[i] = (double) p[k + 3];
			P4[i] = (double) p[k + 4];
			P5[i] = (double) p[k + 5];
			P6[i] = (double) p[k + 6];
			P7[i] = (double) p[k + 7];
			P8[i] = (double) p[k + 8];
			P9[i] = (double) p[k + 9];
			P10[i] = (double) p[k + 10];
			P11[i] = (double) p[k + 11];
			P12[i] = (double) p[k + 12];
			P13[i] = (double) p[k + 13];
			P14[i] = (double) p[k + 14];
		}
		for (int i = L; i < Lpad; i++) {
			P0[i] = 0; P1[i] = 0; P2[i] = 0; P3[i] = 0; P4[i] = 0; P5[i] = 0; P6[i] = 0; P7[i] = 0; P8[i] = 0; P9[i] = 0; P10[i] = 0; P11[i] = 0; P12[i] = 0; P13[i] = 0; P14[i] = 0;
		}
		for (int j = 0; j < N; j++) cr[j] = (double) coalescent_rates[j];
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
		{ double s = 0; for (int i = 0; i < L; i++) s += P13[i]; sS[13] = s; }
		{ double s = 0; for (int i = 0; i < L; i++) s += P14[i]; sS[14] = s; }

		while (duration > 0) {
			iterations++;
			pass1(Lpad);
			pass2(Lpad);
			pass2b(Lpad);
			double timeStep = FastMath.min(FastMath.cbrt(epsilon * 6 / maxDDD), FastMath.min(duration, max_step));
			double timeStepSquare = timeStep * timeStep * 0.5;
			if (needsSearch(Lpad, (double) timeStep, (double) timeStepSquare)) {
				// same sequential search as Euler2ndOrder.updateP, in the original element order
				final double[][] PA = {P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14}, DA = {D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14}, DDA = {DD0, DD1, DD2, DD3, DD4, DD5, DD6, DD7, DD8, DD9, DD10, DD11, DD12, DD13, DD14};
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
								P14[L - 1] = Double.NEGATIVE_INFINITY;
								break search;
							}
						}
					}
				}
			}
			if (P14[L - 1] == Double.NEGATIVE_INFINITY) {
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
			p[k + 5] = P5[i];
			p[k + 6] = P6[i];
			p[k + 7] = P7[i];
			p[k + 8] = P8[i];
			p[k + 9] = P9[i];
			p[k + 10] = P10[i];
			p[k + 11] = P11[i];
			p[k + 12] = P12[i];
			p[k + 13] = P13[i];
			p[k + 14] = P14[i];
		}
		p[last] = logPSlot;
	}
}
