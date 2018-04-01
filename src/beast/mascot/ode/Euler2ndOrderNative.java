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


	
	public static void main(String[] args) {
		loadLibrary();
		Euler2ndOrderNative x = new Euler2ndOrderNative();
		x.setup(123, 2, 0.0001, Double.MAX_VALUE);
	}
}
