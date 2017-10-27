package beast.mascot.operators;

import java.text.DecimalFormat;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;




@Description("A random walk operator that selects a random dimension of the real parameter and perturbs the value a " +
        "random amount within +/- windowSize.")
public class RealRandomWalkOperator extends Operator {
    final public Input<Double> windowSizeOnInput =
            new Input<>("windowSizeOn", "the size of the window both up and down when using uniform interval OR standard deviation when using Gaussian", Input.Validate.REQUIRED);
    final public Input<Double> windowSizeOffInput =
            new Input<>("windowSizeOff", "the size of the window both up and down when using uniform interval OR standard deviation when using Gaussian", Input.Validate.REQUIRED);
    final public Input<RealParameter> parameterInput =
            new Input<>("parameter", "the parameter to operate a random walk on.", Validate.REQUIRED);
    final public Input<Boolean> useGaussianInput =
            new Input<>("useGaussian", "Use Gaussian to move instead of uniform interval. Default false.", true);
    
    final public Input<BooleanParameter> indicatorInput =
            new Input<>("indicator", "Defines which parameters to scale");
   
    
    double windowSizeOn = 1;
    double windowSizeOff = 1;
    boolean useGaussian;

    @Override
	public void initAndValidate() {
        windowSizeOn = windowSizeOnInput.get();
        windowSizeOff = windowSizeOffInput.get();
        useGaussian = useGaussianInput.get();
    }

    /**
     * override this for proposals,
     * returns log of hastingRatio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        RealParameter param = parameterInput.get(this);
        for (int i = 0; i < param.getDimension(); i++){
	        double value = param.getValue(i);
	        double newValue = value;
	        double windowSize = 0.0;
	        
        	if (indicatorInput.get().getArrayValue(i) > 0.5)
        		windowSize = windowSizeOn;
        	else
        		windowSize = windowSizeOff;
        		

	        if (useGaussian) {
	            newValue += Randomizer.nextGaussian() * windowSize;
	        } else {
	            newValue += Randomizer.nextDouble() * 2 * windowSize - windowSize;
	        }
	
	        if (newValue < param.getLower() || newValue > param.getUpper()) {
	            return Double.NEGATIVE_INFINITY;
	        }
	        if (newValue == value) {
	            // this saves calculating the posterior
	            return Double.NEGATIVE_INFINITY;
	        }
	
	        param.setValue(i, newValue);
        }
        return 0.0;
    }


    @Override
    public double getCoercableParameterValue() {
        return windowSizeOn;
    }

    @Override
    public void setCoercableParameterValue(double value) {
        windowSizeOn = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
//        double delta = calcDelta(logAlpha);
//
//        delta += Math.log(windowSize);
//        windowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = windowSizeOn * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }
} // class RealRandomWalkOperator