package beast.mascot.ode;

import beast.core.util.Log;

/** A native implementation of Euler2ndOrder for Mascot **/
public class Euler2ndOrderNative implements Euler2ndOrderBase {

	public static boolean loadLibrary() {
		try {
			System.loadLibrary("mascot");
			Log.warning("Loaded mascot native library");
			return true;
		} catch (java.lang.UnsatisfiedLinkError e) {
			Log.warning("Mascot native library not loaded");
			return false;
		}
	}

	@Override
	native public void setup(int maxSize, int states, double epsilon,
			double max_step);
	
	@Override
	native public void init(double[] migration_rates, double[] coalescent_rates, int lineages);

	@Override
	native public void initWithIndicators(double[] migration_rates, int[] indicators, double[] coalescent_rates,
			int lineages);
	
	@Override
	native public void calculateValues(double duration, double[] p, int length);
	
	@Override
	native public void initAndcalculateValues(double[] migration_rates, double[] coalescent_rates,
			int lineages, double duration, double[] p, int length);

	@Override
	native public void initAndcalculateValues(int ratesInterval, int lineages, double duration, double[] p, int length);
	
	@Override
	public void setUpDynamics(double[][] coalescentRates, double[][] migrationRates, int[][] indicators,
			double[] nextRateShift) {
		int n = coalescentRates.length;
		int n2 = coalescentRates[0].length;
		int n3 = migrationRates[0].length;
		double [] coalescentRatesArray = new double[n*n2];
		double [] migrationRatesArray = new double[n*n3];
		for (int i = 0; i < n; i++) {
			System.arraycopy(coalescentRates[i], 0, coalescentRatesArray, i*n2, n2);
		}
		for (int i = 0; i < n; i++) {
			System.arraycopy(migrationRates[i], 0, migrationRatesArray, i*n3, n3);
		}
		if (indicators[0] == null) {
			setUpDynamics(coalescentRatesArray, migrationRatesArray, nextRateShift);
		} else {
			int n4 = indicators[0].length;
			int [] indicatorsArray = new int[n*n4];
			for (int i = 0; i < n; i++) {
				System.arraycopy(indicators[i], 0, indicatorsArray, i*n4, n4);
			}
			setUpDynamics(coalescentRatesArray, migrationRatesArray, indicatorsArray, nextRateShift);
		}
	}

	native private void setUpDynamics(double[] coalescentRates, double[] migrationRates, double[] nextRateShift);
	native private void setUpDynamics(double[] coalescentRates, double[] migrationRates, int[] indicators,
			double[] nextRateShift);

	
	public static void main(String[] args) {
		loadLibrary();
		Euler2ndOrderNative x = new Euler2ndOrderNative();
		x.setup(123, 2, 0.0001, Double.MAX_VALUE);
	}
}
