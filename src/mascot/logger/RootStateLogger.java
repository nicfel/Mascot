package mascot.logger;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import mascot.distribution.Mascot;

import java.io.PrintStream;

/**
 * @author Nicola Felix Mueller
 */
@Description("logs the state of the root, i.e. the probability the root being in "+ 
			"any of m states based on the MASCOT density class")
public class RootStateLogger extends CalculationNode implements Loggable {
	public Input<Mascot> mascotInput = new Input<Mascot>(
			"mascot",
			"");

	
	private int states;
	@Override
	public void init(PrintStream out) {
		double[] RootStates = mascotInput.get().getRootState();
		states = RootStates.length;
		for (int i = 0 ; i < states; i++){
			out.print("RootProbability." + i + "\t");
		}

	}

	@Override
	public void log(long sample, PrintStream out) {
		double[] RootStates = mascotInput.get().getRootState();
		states = RootStates.length;
		for (int i = 0 ; i < states; i++){
			out.print(RootStates[i] + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {	
	}

	@Override
	public void initAndValidate() {	
	}	
	

	
}
