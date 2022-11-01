package beast.mascot.parameterDynamics;

import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.mascot.dynamics.RateShifts;
import beast.mascot.parameterDynamics.EffectivePopulationSizeDynamics;


/**
 * @author Nicola F. Mueller
 */
@Description("Populaiton function with values at certain time points that are interpolated in between. Parameter has to be in log space")
public class Skygrowth extends EffectivePopulationSizeDynamics {
	
    final public Input<RealParameter> NeInput = new Input<>("logNe",
            "Nes over time in log space", Input.Validate.REQUIRED);
    final public Input<RateShifts> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);

    //
    // Public stuff
    //
    RealParameter Ne;
    RateShifts rateShifts;
    
    boolean NesKnown = false;
    double[] growth;
    double[] growth_stored;


    @Override
	public void initAndValidate() {
    	Ne = NeInput.get();    	    	
    	rateShifts = rateShiftsInput.get();
    	Ne.setDimension(rateShifts.getDimension()+1);
    	growth = new double[rateShifts.getDimension()];
    	recalculateNe();
		isTime = true;

    }


	@Override
	public double getNeTime(double t) {	
//		if (!NesKnown)
//			recalculateNe();

		int intervalnr = getIntervalNr(t);
		if (intervalnr>=rateShifts.getDimension()) {
			return Double.POSITIVE_INFINITY;
		}
		double timediff = t;
		if (intervalnr>0)
			timediff -= rateShifts.getValue(intervalnr-1);
				
		return Math.exp(Ne.getArrayValue(intervalnr)-growth[intervalnr]*timediff);
	}


	
	private int getIntervalNr(double t) {
		// check which interval t + offset is in
		for (int i = 0; i < rateShifts.getDimension(); i++)
			if (t<rateShifts.getValue(i))
				return i;
		
		// after the last interval, just keep using the last element
		return rateShifts.getDimension();					
	}

	
	// computes the Ne's at the break points
	private void recalculateNe() {
		growth = new double[rateShifts.getDimension()];
		double curr_time = 0.0;
		for (int i = 1; i < Ne.getDimension(); i++) {
			growth[i-1] = (Ne.getArrayValue(i-1)- Ne.getArrayValue(i))/(rateShifts.getValue(i-1)-curr_time);
			curr_time = rateShifts.getValue(i-1);
		}
		NesKnown = true;
	}

	@Override
	public boolean requiresRecalculation() {
		recalculateNe();
		return super.requiresRecalculation();
	}
	
	@Override
	public void store() {
		growth_stored = new double[growth.length];
		System.arraycopy(growth, 0, growth_stored, 0, growth.length);
		super.store();
	}
	
	@Override
	public void restore() {
		System.arraycopy(growth_stored, 0, growth, 0, growth_stored.length);
		super.restore();
	}


	@Override
	public void recalculate() {
		// TODO Auto-generated method stub
		
	}


	@Override
	public boolean isDirty() {
		if (Ne.isDirty(0))
			return true;
		
		return false;
	}


	
}