package beast.mascot.dynamics;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.mascot.distribution.StructuredTreeIntervals;

@Description("Calculates coalescent rates based on SIS models")
public class SIS extends Dynamics {
	
    public Input<RealParameter> SInput = new Input<>("S", "input of effective population sizes", Validate.REQUIRED);
    public Input<RealParameter> initialInput = new Input<>("initialIntroducation", "input of effective population sizes", Validate.REQUIRED);
    public Input<RealParameter> infectionInput = new Input<>("infection", "input of effective population sizes", Validate.REQUIRED);
    public Input<RealParameter> recoveryInput = new Input<>("recovery", "input of effective population sizes", Validate.REQUIRED);
    public Input<RealParameter> forwardMigrationInput = new Input<>("forwardMigration", "input of effective population sizes", Validate.REQUIRED);    
    public Input<StructuredTreeIntervals> treeIntervalsInput = new Input<StructuredTreeIntervals>("structuredTreeIntervals", "Structured Intervals for a phylogenetic beast tree", Validate.REQUIRED);
	
    // used to temporarly store the calculated mean coalescent rates in an interval
    private double[] S;
    private double[] I;
    private double[] timepoints;
    private int[] type;
	
	@Override
	public void initAndValidate() {		
	}	
	
	@Override
	public void recalculate() {
		
		HashMap<Integer, Double> introductions = new HashMap<>();
		for (int i = 0; i < initialInput.get().getDimension(); i++)
			introductions.put(i, initialInput.get().getArrayValue(i));
		
		
		// get the order of introduction
		Double[] initialIntroductions = initialInput.get().getValues();
		
		ArrayList<Double> intros = new ArrayList<Double>();
		for(double d : initialIntroductions) intros.add(d);

		// sort the introductions
		Collections.sort((List<T>) introductions);		
		for (int i = intros.size()-1; i > 0; i--)
			intros.set(i, intros.get(i) - intros.get(i-1));
				
		int introNumber = 0, interval = 0, t = 0;

		// get the next intervals
		timepoints = new double[treeIntervalsInput.get().getIntervalCount()+dimensionInput.get()];				
		type = new int[treeIntervalsInput.get().getIntervalCount()+dimensionInput.get()];				
		
		double nextIntroduction = intros.get(introNumber);
		double nextInterval = treeIntervalsInput.get().getInterval(interval);
		
		while (t < timepoints.length){
			if (nextInterval < nextIntroduction){
				interval++;
				timepoints[t] = nextInterval;
				type[t] = Integer.MAX_VALUE;
				t++;
				nextIntroduction -= nextInterval;
				try{
					nextInterval = treeIntervalsInput.get().getInterval(interval);
				}catch (Exception e){
					break;
				}
			}else{
				introNumber++;
				timepoints[t] = nextIntroduction;
				type[t] = 
				t++;
				nextInterval -= nextIntroduction;
				nextIntroduction = intros.get(introNumber);				
			}			
		}
		
		System.out.println(Arrays.toString(timepoints));
		System.exit(0);
		setDynamicsKnown();
	}
	
	@Override
	public double getInterval(int i) {
		return treeIntervalsInput.get().getInterval(i);
	}

	@Override
	public boolean intervalIsDirty(int i) {
		boolean intervalIsDirty = false;
		for (int j = 0; j < dimensionInput.get(); j++)
			if (SInput.get().isDirty(j))
				intervalIsDirty = true;
		for (int j = 0; j < dimensionInput.get(); j++)
			if (initialInput.get().isDirty(j))
				intervalIsDirty = true;
		for (int j = 0; j < dimensionInput.get(); j++)
			if (infectionInput.get().isDirty(j))
				intervalIsDirty = true;
		for (int j = 0; j < dimensionInput.get(); j++)
			if (recoveryInput.get().isDirty(j))
				intervalIsDirty = true;
		for (int j = 0; j < forwardMigrationInput.get().getDimension(); j++)
			if (forwardMigrationInput.get().isDirty(j))
				intervalIsDirty = true;


		return intervalIsDirty;
	}

	@Override
	public double[][] getBackwardsMigration(int interval) {
		if(!areDynamicsKnown())
			recalculate();
		
		double[][] migrationRate = new double[dimensionInput.get()][dimensionInput.get()];
//		
//		for (int a = 0; a < dimensionInput.get(); a++){
//			for (int b = 0; b < dimensionInput.get(); b++){
//				if (a!=b){
//					migrationRate[a][b] 
//							= forwardMigrationInput.get().getArrayValue(a*(dimensionInput.get()-1) + b) 
//								* coalRate[a]/coalRate[b];
//				}	
//			}
//		}
//
		return null;
	}

	@Override
	public double[] getCoalescentRate(int interval) {
		if(!areDynamicsKnown())
			recalculate();
		
		double[] coalRate = new double[dimensionInput.get()];
//		double currentTime = 0.0;
//		
//		if (interval>0)
//			currentTime = treeIntervalsInput.get().getInterval(interval-1);
//		
//		double nextTime = treeIntervalsInput.get().getInterval(interval);
//		
//		// calculate the mean Ne's for each interval
//		for (int i = 0; i < dimensionInput.get(); i++){
//			coalRate[i] = 1/(presentDayNeInput.get().getArrayValue(i)
//					/((currentTime - nextTime)
//							*growthRatesInput.get().getArrayValue(i))
//								*(Math.exp(-growthRatesInput.get().getArrayValue(i)*currentTime) 
//										- Math.exp(-growthRatesInput.get().getArrayValue(i)*nextTime)));
//		}
		return coalRate;
	}
	

}
