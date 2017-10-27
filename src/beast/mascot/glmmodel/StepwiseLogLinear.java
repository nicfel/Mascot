package beast.mascot.glmmodel;

import java.io.PrintStream;

public class StepwiseLogLinear extends GLMmodel {


	@Override
	public void initAndValidate() {
		// set the dimension of the scalers, indicators and potentially the error term
    	scalerInput.get().setDimension(covariatesInput.get().size());
    	indicatorInput.get().setDimension(covariatesInput.get().size());
    	if (errorInput.get()!=null)
    		errorInput.get().setDimension(covariatesInput.get().get(0).getDimension());

		// ensure that all the entries in all the log covariates sum to 0
    	for (int i = 0; i < covariatesInput.get().size(); i++){
    		double paramsum = 0;
    		for (int j = 0; j < covariatesInput.get().get(i).getDimension(); j++)
    			paramsum += Math.log(covariatesInput.get().get(i).getArrayValue(j));
    			
    			paramsum /= covariatesInput.get().get(i).getDimension();
    		for (int j = 0; j < covariatesInput.get().get(i).getDimension(); j++)
    			 covariatesInput.get().get(i).setValue(j, 
    					 Math.log(covariatesInput.get().get(i).getArrayValue(j)) - paramsum);
    	}

	}

	@Override
	public double[] getRates() {
    	double[] logrates = new double[covariatesInput.get().get(0).getDimension()];
    	
    	for (int j = 0; j < covariatesInput.get().get(0).getDimension(); j++)
    		logrates[j] = 0;
    	    	
		for (int j = 0; j < covariatesInput.get().size(); j++){
			if (indicatorInput.get().getArrayValue(j) > 0.0)
				for (int k = 0; k < covariatesInput.get().get(j).getDimension(); k++)
					logrates[k] += scalerInput.get().getArrayValue(j)
						*covariatesInput.get().get(j).getArrayValue(k);
		}
		
    	if (errorInput.get()!=null)
    		for (int k = 0; k < errorInput.get().getDimension(); k++)
    			logrates[k] += errorInput.get().getArrayValue(k);
    	// normalize
    	double ratessum = 0;
    	for (int i = 0; i < logrates.length; i++)
    		ratessum += logrates[i];
    	
    	for (int i = 0; i < logrates.length; i++)
    		logrates[i] -= ratessum/logrates.length ;
    	
    	double[] rates = new double[covariatesInput.get().get(0).getDimension()];
   	
		for (int k = 0; k < covariatesInput.get().get(0).getDimension(); k++)
			rates[k] = clockInput.get().getArrayValue()*Math.exp(logrates[k]);

    	return rates;
	}

	@Override
	public void init(PrintStream out) {
		for (int i = 0 ; i < scalerInput.get().getDimension(); i++){
			out.print(String.format("%sscaler.%s\t", getID(), covariatesInput.get().get(i).getID()));
		}
		
	}

	@Override
	public void log(int sample, PrintStream out) {
		for (int i = 0 ; i < scalerInput.get().getDimension(); i++){
			if (indicatorInput.get().getArrayValue(i) > 0.5)
				out.print(scalerInput.get().getArrayValue(i) +"\t");
			else
				out.print("0.0\t");
		}
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

}
