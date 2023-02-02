package mascot.parameterdynamics;

import java.util.List;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import mascot.dynamics.RateShifts;


/**
 * @author Nicola F. Mueller
 */
@Description("Populaiton function with values at certain time points that are interpolated in between. Parameter has to be in log space")
public class NeSplineInterpolation extends NeDynamics {
	
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
    PolynomialSplineFunction psf;
    SplineInterpolator si;


    @Override
	public void initAndValidate() {
    	Ne = NeInput.get();    	    	
    	rateShifts = rateShiftsInput.get();
    	Ne.setDimension(rateShifts.getDimension()+1);
    	recalculateNe();
		isTime = true;
		si = new SplineInterpolator();
	

    }


	@Override
	public double getNeTime(double t) {	
//		if (!NesKnown)
//			recalculateNe();

		int intervalnr = getIntervalNr(t);
		if (intervalnr>=rateShifts.getDimension()) {
			return Math.exp(Ne.getArrayValue(Ne.getDimension()-1));
		}
		double timediff = t;
		if (intervalnr>0)
			timediff -= rateShifts.getValue(intervalnr-1);
				
		return 0;
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
		double[] nes = Ne.getDoubleValues();
		double[] rates = new double[nes.length];
		for (int i = 0; i < rates.length; i++) {
			rates[i] = rateShifts.rateShifts[i];
		}

		psf = si.interpolate(rates, nes);
		NesKnown = true;
	}

	@Override
	public boolean requiresRecalculation() {
		recalculateNe();
		return super.requiresRecalculation();
	}
	
	@Override
	public void store() {
//		growth_stored = new double[growth.length];
//		System.arraycopy(growth, 0, growth_stored, 0, growth.length);
		super.store();
	}
	
	@Override
	public void restore() {
//		System.arraycopy(growth_stored, 0, growth, 0, growth_stored.length);
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