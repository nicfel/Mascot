package beast.mascot.glmmodel;

import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.Parameter.Base;
import beast.core.parameter.RealParameter;

public class CovariateList extends BEASTObject {
    public Input<List<Covariate>> covariatesInput = new Input<>("covariates", "input of covariates", new ArrayList<>(), Validate.REQUIRED);

    
    List<Covariate> covariates;
    
    public CovariateList() {
    	covariates = covariatesInput.get();
    }
   
	@Override
	public void initAndValidate() {
    	covariates = covariatesInput.get();
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
	

}
