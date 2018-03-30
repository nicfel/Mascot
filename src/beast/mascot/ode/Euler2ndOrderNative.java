package beast.mascot.ode;


/** A native implementation of Euler2ndOrder for Mascot **/
public class Euler2ndOrderNative implements Euler2ndOrderBase {

	public static void loadLibrary() {
		System.loadLibrary("mascot");
	}

	@Override
	native public void setup(int maxSize);
	
	@Override
	native public void init(double[] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon,
			double max_step);

	@Override
	native public void initWithIndicators(double[] migration_rates, int[] indicators, double[] coalescent_rates,
			int lineages, int states, double epsilon, double max_step);
	
	@Override
	native public void calculateValues(double duration, double[] p, int length);


	
	public static void main(String[] args) {
		loadLibrary();
		Euler2ndOrderNative x = new Euler2ndOrderNative();
		x.setup(123);
	}
}
