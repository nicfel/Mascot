package mascot.util;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;


@Description("calculates the differences between the entries of a vector")
public class Final extends CalculationNode implements Function {
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
    	return functionInput.get().getArrayValue(functionInput.get().getDimension()-1);
     }


    @Override
    public double getArrayValue(int dim) {
    	return functionInput.get().getArrayValue(functionInput.get().getDimension()-1);
    }

} // class Sum
