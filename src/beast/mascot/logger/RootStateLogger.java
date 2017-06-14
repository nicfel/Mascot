package beast.mascot.logger;

import java.io.PrintStream;

import org.jblas.DoubleMatrix;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.mascot.distribution.Mascot;

/**
 * @author Nicola Felix Mueller
 */
@Description("logs the state of the root, i.e. the probability the root being in "+ 
			"any of m states based on the MASCOT density class")
public class RootStateLogger extends CalculationNode implements Loggable {
	public Input<Mascot> mascotInput = new Input<Mascot>(
			"mascot",
			"A deterministic epidemiological model");

	
	private int states;
	@Override
	public void init(PrintStream out) {
		DoubleMatrix RootStates = mascotInput.get().getRootState();
		states = RootStates.length;
		for (int i = 0 ; i < states; i++){
			out.print("RootProbability." + i + "\t");
		}

	}

	@Override
	public void log(int sample, PrintStream out) {
		DoubleMatrix RootStates = mascotInput.get().getRootState();
		states = RootStates.length;
		for (int i = 0 ; i < states; i++){
			out.print(RootStates.get(i) + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {	
	}

	@Override
	public void initAndValidate() {	
	}	
	

	
}
