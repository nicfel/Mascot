package mascot.util;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;


@Description("calculates the differences between the entries of a vector")
public class Difference extends CalculationNode implements Function {
    final public Input<Function> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double[] difference;
    double[] storedDifference;

    @Override
    public void initAndValidate() {
    	difference = new double[functionInput.get().getDimension()];
    	storedDifference = new double[functionInput.get().getDimension()];
    }

    @Override
    public int getDimension() {
        return difference.length;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return difference[0];
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
        for (int i = 1; i < functionInput.get().getDimension(); i++) {
        	difference[i-1] = functionInput.get().getArrayValue(i-1)-functionInput.get().getArrayValue(i);
        }
        
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return difference[dim];
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
    	System.arraycopy(difference, 0, storedDifference, 0, difference.length);
        super.store();
    }

    @Override
    public void restore() {
    	double [] tmp = storedDifference;
    	storedDifference = difference;
    	difference = tmp;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }
} // class Sum
