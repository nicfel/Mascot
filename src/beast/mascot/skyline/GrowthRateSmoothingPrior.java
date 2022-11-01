package beast.mascot.skyline;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.mascot.dynamics.RateShifts;
import beast.math.distributions.ParametricDistribution;

public class GrowthRateSmoothingPrior extends Distribution {
	
    public Input<RealParameter> NeLogInput = new Input<>(
    		"NeLog", "input of effective population sizes");        
    
    final public Input<ParametricDistribution> distInput = new Input<>("distr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.", 
    		Validate.REQUIRED);
    
    final public Input<ParametricDistribution> initDistrInput = new Input<>("initialDistr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.",
    		Input.Validate.OPTIONAL);
        
    final public Input<ParametricDistribution> finalDistrInput = new Input<>("finalDistr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.", 
    		Input.Validate.OPTIONAL);
    
    final public Input<ParametricDistribution> meanDistrInput = new Input<>("meanDistr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.", 
    		Input.Validate.OPTIONAL);
    
    public Input<RateShifts> rateShiftsInput = new Input<>(
    		"rateShifts", "timing of the rate shifts", Validate.REQUIRED);   


   
    
    private RealParameter NeLog;
    
    protected ParametricDistribution dist;
    protected ParametricDistribution initDistr;
    protected ParametricDistribution finalDistr;
    protected ParametricDistribution meanDistr;
    
    RateShifts rateShifts;
    
    @Override
    public void initAndValidate() {
    	NeLog = NeLogInput.get();    	
        dist = distInput.get();
        rateShifts = rateShiftsInput.get();
        
        if (initDistrInput.get()!=null)
        	initDistr = initDistrInput.get();
        if (finalDistrInput.get()!=null)
        	finalDistr = finalDistrInput.get();
        if (meanDistrInput.get()!=null)
        	meanDistr = meanDistrInput.get();
        

    }


	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
    public double calculateLogP() {
        logP = 0;
        
        double[] growthRates = new double[NeLog.getDimension()-1];
               
        
        //loop over all time points
    	for (int j = 1; j < NeLog.getDimension(); j++){
    		double timediff = rateShifts.getValue(j) - rateShifts.getValue(j-1);
    		double logdiff = NeLog.getArrayValue(j) - NeLog.getArrayValue(j-1); 
    		growthRates[j-1] = logdiff/timediff;    		
    	}
    	
    	
    	for (int j = 1; j < growthRates.length; j++){
    		double diff = growthRates[j]-growthRates[j-1];
        	logP += dist.logDensity(diff);
    	}
        
        // add contribution from first or last entry
        if (initDistrInput.get()!=null)
    		logP += initDistr.logDensity(growthRates[0]);
        if (finalDistrInput.get()!=null) {
    		logP += finalDistr.logDensity(growthRates[growthRates.length-1]);
        }
        
        if (meanDistrInput.get()!=null) {
        	double mean=0.0;
        	for (int j = 0; j < growthRates.length; j++){
        		mean += growthRates[j];    		
        	}
        	mean /= growthRates.length;
    		logP += meanDistr.logDensity(mean);
        }
        
        return logP;
    }
}
