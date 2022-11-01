package beast.mascot.dynamics;



import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.mascot.dynamics.Dynamics;
import beast.mascot.dynamics.RateShifts;
import beast.mascot.parameterDynamics.EffectivePopulationSizeDynamics;
import beast.mascot.parameterDynamics.StructuredSkygrid;


@Description("Wrapper class that takes parametric and non parametric dynamics as input.")
public class DynamicEffectivePopulationSizesInterpolate extends Dynamics implements Loggable {	
    
    
    public Input<List<EffectivePopulationSizeDynamics>> parametricFunctionInput = new Input<>(
    		"NeDynamics", "input of the log effective population sizes", new ArrayList<>());    
    
    public Input<RealParameter> b_mInput = new Input<>(
    		"backwardsMigration", "input of backwards in time migration rates");    
    public Input<RealParameter> f_mInput = new Input<>(
    		"forwardsMigration", "input of backwards in time migration rates", Validate.XOR, b_mInput);  
    
    public Input<RateShifts> rateShiftsInput = new Input<>(
    		"rateShifts", "timing of the rate shifts", Validate.REQUIRED);   
    
	public Input<Double> maxRateInput = new Input<>(
			"maxRate", "maximum rate used for integration", Double.POSITIVE_INFINITY);
	
	public Input<Integer> nrInterpolationPointsInput = new Input<>(
			"nrInterpolationPoints", "maximum rate used for integration", 5);

	double[] intTimes;	

	boolean isForward = false;
	
	RealParameter migration;
	
	List<EffectivePopulationSizeDynamics> parametricFunction;
	
	
	int nrPoints;	
    @Override
    public void initAndValidate() {
    	
    	super.initAndValidate();
    	
    	if (dimensionInput.get()<1)
    		dimensionInput.set(getNrTypes());
    	
    	if (f_mInput.get()!=null){
    		migration = f_mInput.get();
    		isForward = true;
    	}else{
    		migration = b_mInput.get();
    	}
    	// check the dimension of the migration rates
    	if (dimensionInput.get()*(dimensionInput.get()-1)!=migration.getDimension() || fromBeautiInput.get()){
    		System.err.println("the dimension of " + migration.getID() + " is set to " + (dimensionInput.get()*(dimensionInput.get()-1)));
    		migration.setDimension(dimensionInput.get()*(dimensionInput.get()-1));
    	}
    	
    	parametricFunction = parametricFunctionInput.get();
    	    	
    	// check that all the inputs have the correct dimension
    	if (dimensionInput.get()!=parametricFunction.size()){
    		throw new IllegalArgumentException("the number of parametric functions given as Input does not match the dimension/number of states");
    	}
    	
		// initialize the intervals
		nrPoints = nrInterpolationPointsInput.get()+1;
		intTimes = new double[rateShiftsInput.get().getDimension()*nrPoints];
		intTimes[0] = rateShiftsInput.get().getValue(0);				
		for (int i = 1; i < rateShiftsInput.get().getDimension(); i++) {
			for (int j = 0; j < nrPoints; j++) {
				intTimes[i] = (rateShiftsInput.get().getValue(i)-rateShiftsInput.get().getValue(i-1))/(nrPoints);
			}
		}
		// set the Ne dimension for all skygrid dynamics.
		for (int i = 0; i < parametricFunction.size(); i++)
			parametricFunction.get(i).setNrIntervals(rateShiftsInput.get().getDimension());
    }

    /**
     * Returns the time to the next interval.
     */
    @Override
    public double getInterval(int i) {
    	if (i >= intTimes.length){
     		return Double.POSITIVE_INFINITY;
     	}else{
			return intTimes[i];
     	}
    }   

    @Override
    public double[] getIntervals() {
    	return intTimes;
    }
    
	@Override
    public double[] getCoalescentRate(int i){
		int intervalNr;
    	if (i > intTimes.length)
    		intervalNr = rateShiftsInput.get().getDimension()-2;
    	else
    		intervalNr = i/nrPoints;
    	
    	int rest = i % nrPoints;
    	System.out.println(i + " " + nrPoints + " " + intervalNr + " " + rest);

    	// get the current time as the midpoint of the current intervals
    	double t = rateShiftsInput.get().getIntervalMidpoint(intervalNr);   

		double[] coal = new double[getDimension()];
		for (int j = 0; j < coal.length; j++){
			if (parametricFunction.get(j).isTime)
				coal[j] = Math.min(maxRateInput.get(), 1/parametricFunction.get(j).getNeTime(t));
			else
				coal[j] = Math.min(maxRateInput.get(), 1/parametricFunction.get(j).getNeInterval(intervalNr));
		}
		return coal;
    }
    
	@Override    
    public double[] getBackwardsMigration(int i){
		int intervalNr;
    	if (i > intTimes.length)
    		intervalNr = rateShiftsInput.get().getDimension()-2;
    	else
    		intervalNr = i/nrPoints;
    	
    	int rest = i % nrPoints;
//    	System.out.println(i + " " + nrPoints + " " + intervalNr + " " + rest);

		
    	int dim = dimensionInput.get();
    	int start = intervalNr*dim*(dim-1);
    	
		double[] m = new double[dim * dim];
		
		if (isForward){
			double[] NeInt = new double[dim];
	    	double t = rateShiftsInput.get().getIntervalMidpoint(intervalNr);   	
			
			for (int j = 0; j < NeInt.length; j++){
				if (parametricFunction.get(j).isTime)
					NeInt[j] = parametricFunction.get(j).getNeTime(t);
				else
					NeInt[j] = parametricFunction.get(j).getNeInterval(intervalNr);
			}

			int c=0;
			
			for (int a = 0; a < dim; a++){
				for (int b = 0; b < dim; b++){
					if (a!=b){
						m[b * dim + a] = Math.min(maxRateInput.get(), NeInt[a]*migration.getArrayValue(c)/NeInt[b]);
						c++;
					}
				}
			}

			
		}else{
			int c1=0,c2=0;
			start = 0;
			for (int a = 0; a < dim; a++){
				for (int b = 0; b < dim; b++){
					if (a!=b){
						m[c1] = Math.min(maxRateInput.get(), migration.getArrayValue(start+c2));
						c2++;
					}
					c1++;
				}
			}
		}		
		return m;
    }

//	@Override
//	public void recalculate() {
//		// TODO Auto-generated method stub
//		
//	}    
//	
//	@Override
//	public void init(PrintStream out) {
//		for (int j = 0; j < dimensionInput.get(); j++){
//			for (int i = 0; i < rateShifts.getDimension(); i++){
//				out.print(String.format("Ne.%d.%d\t", j,i));
//			}			
//		}
//	}
//
//	@Override
//	public void log(long sample, PrintStream out) {
//		for (int j = 0; j < dimensionInput.get(); j++){
//			for (int i = 0; i < rateShifts.getDimension(); i++){
//				out.print(Math.exp(Ne.getArrayValue(dimensionInput.get()*i+j))*NeMean.getValue() + "\t");
//			}			
//		}
//	}

    @Override
    public int getEpochCount() {
    	return rateShiftsInput.get().getCorrectedDimension();
    }

	@Override
	public boolean intervalIsDirty(int j) {
		
		for (int i = 0; i < parametricFunction.size(); i++)
			if(parametricFunction.get(i).isDirty())
				return true;
				
		for (int i = 0; i < migration.getDimension(); i++)
			if(migration.isDirty(i))
				return true;

		if(rateShiftsInput.get().somethingIsDirty())
			return true;

		return false;
	}

	@Override
	public void init(PrintStream out) {
		if (isForward){
			for (int a = 0; a < dimensionInput.get(); a++) {
				for (int b = 0; b < intTimes.length; b++) {
					double t = rateShiftsInput.get().getIntervalMidpoint(b);
					double ne;
					if (parametricFunction.get(a).isTime)
						ne = parametricFunction.get(a).getNeTime(t);
					else
						ne = parametricFunction.get(a).getNeInterval(b);

					out.print(String.format("logNe_%s.%d\t", getStringStateValue(a), b));
				}
			}			

			for (int a = 0; a < dimensionInput.get(); a++)
				for (int b = 0; b < dimensionInput.get(); b++)
					if (a!=b)
						out.print(String.format("f_%s.%s_to_%s\t", f_mInput.get().getID(), getStringStateValue(a), getStringStateValue(b)));
			}			
	}
	
	@Override
	public void log(long sample, PrintStream out) {
		int c=0;
		if (isForward){
			for (int a = 0; a < dimensionInput.get(); a++) {
				for (int b = 0; b < intTimes.length; b++) {
					double t = rateShiftsInput.get().getIntervalMidpoint(b);
					double ne;
					if (parametricFunction.get(a).isTime)
						ne = parametricFunction.get(a).getNeTime(t);
					else
						ne = parametricFunction.get(a).getNeInterval(b);

					out.print(Math.log(ne) + "\t");
				}
			}

			for (int a = 0; a < dimensionInput.get(); a++) {
				for (int b = 0; b < dimensionInput.get(); b++) { 
					if (a!=b) {
						out.print(f_mInput.get().getArrayValue(c) + "\t");
						c++;
					}
				}	
			}
		}			
		
	}


	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void recalculate() {
		// TODO Auto-generated method stub
		
	}
	
}