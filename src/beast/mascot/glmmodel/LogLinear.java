package beast.mascot.glmmodel;

import java.io.PrintStream;

import beast.core.Input;

public class LogLinear extends GlmModel {

	@Override
	public void initAndValidate() {
		// set the dimension of the scalers, indicators and potentially the error term
    	scalerInput.get().setDimension(covariateListInput.get().size());
    	indicatorInput.get().setDimension(covariateListInput.get().size());
    	
    	if (errorInput.get()!=null)
    		errorInput.get().setDimension(covariateListInput.get().get(0).getDimension());
    	
    	if (constantErrorInput.get()!=null)
    		if (constantErrorInput.get().getDimension()<1)
    			constantErrorInput.get().setDimension(verticalEntries);
	}
	

	@Override
	public double[] getRates(int i) {
    	double[] logrates = new double[verticalEntries];
    	
    	for (int j = 0; j < logrates.length; j++)
    		logrates[j] = 0;
    	    	
		for (int j = 0; j < covariateListInput.get().size(); j++){
			if (indicatorInput.get().getArrayValue(j) > 0.0){
				for (int k = 0; k < logrates.length; k++){
					logrates[k] += scalerInput.get().getArrayValue(j)
						*covariateListInput.get().get(j).getArrayValue(verticalEntries*i + k);
				}
			}
		}
		
    	if (errorInput.get()!=null)
    		for (int k = 0; k < logrates.length; k++)
    			logrates[k] += errorInput.get().getArrayValue(verticalEntries*i + k);
    	
    	if (constantErrorInput.get()!=null)
    		for (int k = 0; k < logrates.length; k++)
    			logrates[k] += constantErrorInput.get().getArrayValue( k);

    	double[] rates = new double[verticalEntries];
   	
		for (int k = 0; k < verticalEntries; k++){
			rates[k] = clockInput.get().getArrayValue()*Math.exp(logrates[k]);				
		}

		return rates;
	}

	@Override
	public void init(PrintStream out) {
		out.print(String.format("%sClock\t", getID()));

		for (int i = 0 ; i < scalerInput.get().getDimension(); i++){
			out.print(String.format("%sscaler.%s\t", getID(), covariateListInput.get().get(i).getID()));
		}
	}

	@Override
	public void log(int sample, PrintStream out) {
		out.print(clockInput.get().getValue() +"\t");

		for (int i = 0 ; i < scalerInput.get().getDimension(); i++){
			if (indicatorInput.get().getArrayValue(i) > 0.5){
				out.print(scalerInput.get().getArrayValue(i) +"\t");
			}else{
				out.print("0.0\t");
			}
		}
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

}
