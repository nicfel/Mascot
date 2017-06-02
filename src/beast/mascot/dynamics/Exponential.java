package beast.mascot.dynamics;

import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.mascot.distribution.StructuredTreeIntervals;

@Description("Calculate Ne's and backwards in time migration rates for an exponential growth model")
public class Exponential extends Dynamics {
	
    public Input<RealParameter> presentDayNeInput = new Input<>("presentNe", "input of effective population sizes", Validate.REQUIRED);    
    public Input<RealParameter> growthRatesInput = new Input<>("growthRate", "input of effective population sizes", Validate.REQUIRED);    
    public Input<RealParameter> forwardMigrationInput = new Input<>("forwardMigration", "input of effective population sizes", Validate.REQUIRED);    
    public Input<StructuredTreeIntervals> treeIntervalsInput = new Input<StructuredTreeIntervals>("structuredTreeIntervals",
    		"Structured Intervals for a phylogenetic beast tree", Validate.REQUIRED);
	
    // used to temporarly store the calculated mean coalescent rates in an interval
    private double[] coalRate;
    private double oldTime, newTime;
	
	@Override
	public void initAndValidate() {		
	}	
	
	@Override
	public void recalculate() {
	}
	
	@Override
	public double getInterval(int i) {
		return treeIntervalsInput.get().getInterval(i);
	}

	@Override
	public boolean intervalIsDirty(int i) {
		boolean intervalIsDirty = false;
		for (int j = 0; j < dimensionInput.get(); j++)
			if (presentDayNeInput.get().isDirty(j))
				intervalIsDirty = true;
		for (int j = 0; j < dimensionInput.get(); j++)
			if (growthRatesInput.get().isDirty(j))
				intervalIsDirty = true;
		for (int j = 0; j < forwardMigrationInput.get().getDimension(); j++)
			if (forwardMigrationInput.get().isDirty(j))
				intervalIsDirty = true;


		return intervalIsDirty;
	}

	@Override
	public double[][] getBackwardsMigration(int interval) {
		double[][] migrationRate = new double[dimensionInput.get()][dimensionInput.get()];
		
		for (int a = 0; a < dimensionInput.get(); a++){
			for (int b = 0; b < dimensionInput.get(); b++){
				if (a!=b){
					migrationRate[a][b] 
							= forwardMigrationInput.get().getArrayValue(a*(dimensionInput.get()-1) + b) 
								* coalRate[a]/coalRate[b];
				}	
			}
		}

		return migrationRate;
	}

	@Override
	public double[] getCoalescentRate(int interval) {
		coalRate = new double[dimensionInput.get()];
		
		double[] times = treeIntervalsInput.get().getIntervalStartEnd(interval);
//		System.out.println("times " + Arrays.toString(times));
		
		// calculate the mean Ne's for each interval
		for (int i = 0; i < dimensionInput.get(); i++){
			coalRate[i] = 1/(presentDayNeInput.get().getArrayValue(i)
					/((times[1] - times[0])
							*growthRatesInput.get().getArrayValue(i))
								*(Math.exp(-growthRatesInput.get().getArrayValue(i)*times[0]) 
										- Math.exp(-growthRatesInput.get().getArrayValue(i)*times[1])));
		}
//		System.out.println("coalRate " + Arrays.toString(coalRate));
		return coalRate;
	}
	

}
