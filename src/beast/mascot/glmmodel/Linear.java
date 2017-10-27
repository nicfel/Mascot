package beast.mascot.glmmodel;

import java.io.PrintStream;

public class Linear extends GLMmodel {


	@Override
	public void initAndValidate() {
		// set the dimension of the scalers, indicators and potentially the error term
    	scalerInput.get().setDimension(covariatesInput.get().size());
    	indicatorInput.get().setDimension(covariatesInput.get().size());
    	if (errorInput.get()!=null)
    		errorInput.get().setDimension(covariatesInput.get().get(0).getDimension());

		// ensure that all the entries in all the covariates have a mean of 1
    	for (int i = 0; i < covariatesInput.get().size(); i++){
    		double paramsum = 0;
    		for (int j = 0; j < covariatesInput.get().get(i).getDimension(); j++)
    			paramsum += covariatesInput.get().get(i).getArrayValue(j);
    		for (int j = 0; j < covariatesInput.get().get(i).getDimension(); j++)
    			 covariatesInput.get().get(i).setValue(j, 
    					 covariatesInput.get().get(i).getArrayValue(j)/paramsum
    					 	* covariatesInput.get().get(i).getDimension());
    	}

	}

	@Override
	public double[] getRates() {
    	double[] rates = new double[covariatesInput.get().get(0).getDimension()];
    	
    	for (int j = 0; j < covariatesInput.get().get(0).getDimension(); j++)
    		rates[j] = 0;
    	
		for (int j = 0; j < covariatesInput.get().size(); j++){
			if (indicatorInput.get().getArrayValue(j) > 0.0)
				for (int k = 0; k < covariatesInput.get().get(j).getDimension(); k++)
					rates[k] += scalerInput.get().getArrayValue(j)
						*covariatesInput.get().get(j).getArrayValue(k);
		}
    	if (errorInput.get()!=null)
    		for (int k = 0; k < errorInput.get().getDimension(); k++)
    			rates[k] += errorInput.get().getArrayValue(k);

		for (int k = 0; k < errorInput.get().getDimension(); k++)
			rates[k] *= clockInput.get().getArrayValue();

    	return rates;
	}

	@Override
	public void init(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void log(int sample, PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

}
