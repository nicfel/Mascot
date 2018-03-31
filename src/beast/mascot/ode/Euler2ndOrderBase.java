package beast.mascot.ode;

public interface Euler2ndOrderBase {

	public void setup(int maxSize);
	
	public void init(double[] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step);
	public void initWithIndicators(double[] migration_rates, int[] indicators, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step);
	public void calculateValues(double duration, double[] p, int length);
	default public void initAndcalculateValues(double[] migration_rates, double[] coalescent_rates, int lineages, int states, double epsilon, double max_step, double duration, double[] p, int length) {
		init(migration_rates, coalescent_rates, lineages, states, epsilon, max_step);
		calculateValues(duration, p, length);
	}

}
