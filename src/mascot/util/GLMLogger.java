package mascot.util;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.mascot.dynamics.GLM;

@Description("Logs Ne's and migration rates from a GLM model.")

public class GLMLogger extends CalculationNode implements Loggable {
	public Input<GLM> glmInput = new Input<>(
			"glm", "input of GLM model", Validate.REQUIRED);
	
	public Input<Boolean> NeTimeVariantInput = new Input<>(
			"NeTimeVariant", "input of GLM model", true);
	public Input<Boolean> migTimeVariantInput = new Input<>(
			"migTimeVariant", "input of GLM model", true);


		
	int nrIntervals,dims;
	GLM glm;
	boolean NeTimeVariant, migTimeVariant;
	
	@Override
	public void initAndValidate() {
		glm = glmInput.get();
		nrIntervals = glm.getEpochCount();
		dims = glm.getDimension();
		
		NeTimeVariant = NeTimeVariantInput.get();
		migTimeVariant = migTimeVariantInput.get();
	}

	
	@Override
	public void init(PrintStream out) {
		
		for (int i = 0 ; i < dims; i++){
			if (NeTimeVariant) {
				for (int j = 0; j < nrIntervals; j++) {
					out.print(String.format("Ne.%s.%d\t", glm.getStringStateValue(i), j));
				}
			}else {
				for (int j = 0; j < nrIntervals; j++) {
					out.print(String.format("Ne.%s\t", glm.getStringStateValue(i)));
				}
			}
		}
		
		for (int a = 0 ; a < dims; a++){
			for (int b = 0 ; b < dims; b++){
				if (a!=b) {
					if (migTimeVariant) {
						for (int j = 0; j < nrIntervals; j++) {
							out.print(String.format("mig.%s_to_%s.%d\t", glm.getStringStateValue(a), glm.getStringStateValue(b), j));
						}
					}else{
						out.print(String.format("mig.%s_to_%s\t", glm.getStringStateValue(a), glm.getStringStateValue(b)));
					}
				}
			}
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0 ; i < dims; i++){
			if (NeTimeVariant) {
				for (int j = 0; j < nrIntervals; j++) {
					out.print(glm.getNe(i,j)  + "\t");
				}
			}else {
				out.print(glm.getNe(i,0)  + "\t");				
			}
		}
		
		for (int a = 0 ; a < dims; a++){
			for (int b = 0 ; b < dims; b++){
				if (a!=b) {
					if (migTimeVariant) {
						for (int j = 0; j < nrIntervals; j++) {
							out.print(glm.getMig(a,b,j)  + "\t");
						}
					}else {
						out.print(glm.getMig(a,b,0)  + "\t");
					}
				}
			}
		}		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}



}
