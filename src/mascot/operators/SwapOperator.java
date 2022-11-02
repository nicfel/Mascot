package mascot.operators;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;



@Description("A generic operator swapping a one or more pairs in a multi-dimensional parameter")
public class SwapOperator extends Operator {
    final public Input<RealParameter> realparameterInput = new Input<>("parameter", "a real parameter to swap individual values for", Validate.OPTIONAL);
    final public Input<Integer> howManyInput = new Input<>("howMany", "number of items to swap, default 1, must be less than half the dimension of the parameter", 1);


    Parameter<?> parameter;

    @Override
    public void initAndValidate() {
    	parameter = realparameterInput.get();
    }

    @Override
    public double proposal() {
        int left = Randomizer.nextInt(parameter.getDimension());            
    	int right = Randomizer.nextInt(parameter.getDimension());
        
    	while (left==right)
    		right = Randomizer.nextInt(parameter.getDimension());
        
        parameter.swap(left, right);
        
        return 0.0;
    }

}
