package beast.mascot.operators;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;



@Description("A generic operator swapping a one or more pairs in a multi-dimensional parameter")
public class BooleanSwapOperator extends Operator {
    final public Input<BooleanParameter> boolparameterInput = new Input<>("indicator", "an indicator parameter to swap individual values for", Validate.REQUIRED);
    final public Input<RealParameter> realparameterInput = new Input<>("parameter", "a real parameter to swap individual values for", Validate.OPTIONAL);
    final public Input<Integer> howManyInput = new Input<>("howMany", "number of items to swap, default 1, must be less than half the dimension of the parameter", 1);


    int howMany;
    Parameter<?> indicator;
    Parameter<?> parameter;
    private List<Integer> masterList = null;

    @Override
    public void initAndValidate() {
    	indicator = boolparameterInput.get();
        if (realparameterInput.get()!=null){
	        parameter = realparameterInput.get();
	        if (indicator.getDimension()!=parameter.getDimension()){
	            throw new IllegalArgumentException("indicator and parameter have different dimensions");
	        }
        }
        
        howMany = howManyInput.get();
        if (howMany * 2 > indicator.getDimension()) {
            throw new IllegalArgumentException("howMany it too large: must be less than half the dimension of the parameter");
        }

        List<Integer> list = new ArrayList<>();
        for (int i = 0; i < indicator.getDimension(); i++) {
            list.add(i);
        }
        masterList = Collections.unmodifiableList(list);
    }

    @Override
    public double proposal() {
        List<Integer> allIndices = new ArrayList<>(masterList);
        int left, right;

        for (int i = 0; i < howMany; i++) {
            left = allIndices.remove(Randomizer.nextInt(allIndices.size()));            
        	right = allIndices.remove(Randomizer.nextInt(allIndices.size()));
            
            // repeat until left and right are different
            if (indicator.getArrayValue(left)==indicator.getArrayValue(right)){
            	return Double.NEGATIVE_INFINITY;
            }            
            
            indicator.swap(left, right);
            if (realparameterInput.get()!=null){
            	parameter.swap(left, right);
            }
        }
        


        return 0.0;
    }

}
