package mascot.util;


import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;


@Description("calculates the differences between the entries of a vector")
public class Mean extends CalculationNode implements Function {
    final public Input<Function> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double mean;
    double storedMean;

    @Override
    public void initAndValidate() {
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return mean;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
    	mean = 0.0;
        for (int i = 0; i < functionInput.get().getDimension(); i++) {
        	mean += functionInput.get().getArrayValue(i);
        }
        
        mean/=functionInput.get().getDimension();
        
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return mean;
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
    	storedMean = mean;
        super.store();
    }

    @Override
    public void restore() {
    	mean = storedMean;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }
} // class Sum
