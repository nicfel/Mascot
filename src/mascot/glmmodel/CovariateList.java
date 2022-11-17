package mascot.glmmodel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.BooleanParameter;

public class CovariateList extends BEASTObject {
    public Input<List<Covariate>> covariatesInput = new Input<>("covariates", "input of covariates", new ArrayList<>(), Validate.REQUIRED);

    public Input<BooleanParameter> transformInput = new Input<>("transform", "whether covariates need transformation");

    List<Covariate> covariates;
    boolean[] timeDependent;

    // needed to read in predictors with BEAUti
    public HashMap<String, Integer> traitToType = new HashMap<>(); 

    public int nrIntervals=1; 

    
    public CovariateList() {
    	covariates = covariatesInput.get();
    }
   
	@Override
	public void initAndValidate() {
    	covariates = covariatesInput.get();
    	if (transformInput.get()!=null)
    		for (int i = 0; i < covariates.size(); i++)
    			if (transformInput.get().getArrayValue(i)>0.5)
    				covariates.get(i).transform();
    			
	}

	public int size() {
		return covariatesInput.get().size();
	}

	public Covariate get(int index) {
		if (covariatesInput.get()==null)
			return null;
		
		return covariatesInput.get().get(index);
	}
	
	public void add(Covariate cov) {
		 covariatesInput.get().add(cov);
	}

	public String initMigrationFromRawValues(int index) {
		 return covariatesInput.get().get(index).initMigrationFromRawValues(traitToType, nrIntervals);
	}

	public String initNeFromRawValues(int index) {
		 return covariatesInput.get().get(index).initNeFromRawValues(traitToType, nrIntervals);
	}

}
