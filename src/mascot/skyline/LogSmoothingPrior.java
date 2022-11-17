package mascot.skyline;

import java.util.List;
import java.util.Random;

import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.distribution.ParametricDistribution;

public class LogSmoothingPrior extends Distribution {
	
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
    
    private RealParameter NeLog;
    
    protected ParametricDistribution dist;
    protected ParametricDistribution initDistr;
    protected ParametricDistribution finalDistr;
    protected ParametricDistribution meanDistr;
    
    @Override
    public void initAndValidate() {
    	NeLog = NeLogInput.get();    	
        dist = distInput.get();
        
        
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
        
               
        
        //loop over all time points
    	for (int j = 1; j < NeLog.getDimension(); j++){
    		double diff = NeLog.getArrayValue(j) - NeLog.getArrayValue(j-1);    		
    		logP += dist.logDensity(diff);
    	}
        
        // add contribution from first or last entry
        if (initDistrInput.get()!=null)
    		logP += initDistr.logDensity(NeLog.getArrayValue(0));
        if (finalDistrInput.get()!=null) {
    		logP += finalDistr.logDensity(NeLog.getArrayValue(NeLog.getDimension()-1));
        }
        
        if (meanDistrInput.get()!=null) {
        	double mean=0.0;
        	for (int j = 0; j < NeLog.getDimension(); j++){
        		mean += NeLog.getArrayValue(j);    		
        	}
        	mean /= NeLog.getDimension();
    		logP += meanDistr.logDensity(mean);
        }
        
        return logP;
    }
}
