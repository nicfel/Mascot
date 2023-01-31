package mascot.parameterdynamics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.BooleanParameter;

public class NeDynamicsList extends CalculationNode {
    public Input<List<NeDynamics>> neDynamicsInput = new Input<>("neDynamics", "input of covariates", new ArrayList<>(), Validate.REQUIRED);

    List<NeDynamics> neDynamics;

    // needed to read in predictors with BEAUti
    public HashMap<String, Integer> traitToType = new HashMap<>(); 

    public int nrIntervals=1; 

    
    public NeDynamicsList() {
    	neDynamics = neDynamicsInput.get();
    }
   
	@Override
	public void initAndValidate() {
		neDynamics = neDynamicsInput.get();
	}

	public int size() {
		return neDynamicsInput.get().size();
	}

	public NeDynamics get(int index) {
		if (neDynamicsInput.get()==null)
			return null;
		
		return neDynamicsInput.get().get(index);
	}
	
	public void add(NeDynamics dyn) {
		neDynamicsInput.get().add(dyn);
	}

	public void clear() {
		neDynamicsInput.get().clear();
		
	}
}
