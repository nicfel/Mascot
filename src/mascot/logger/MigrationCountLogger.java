package mascot.logger;

import java.io.PrintStream;

import beast.base.inference.CalculationNode;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import cern.colt.Arrays;
import mascot.distribution.MappedMascot;

public class MigrationCountLogger extends CalculationNode implements Loggable {
	public Input<MappedMascot> mappedInput = new Input<MappedMascot>(
			"mappedMascot", "mappedMascot", Input.Validate.REQUIRED);
	
	MappedMascot mapped;
	int[][] counts;
	
	
	@Override
	public void initAndValidate() {
		mapped = mappedInput.get();	
		
	}

	@Override
	public void init(PrintStream out) {
				

		for (int a=0; a < mapped.dynamics.getDimension(); a++) {
			for (int b = 0; b < mapped.dynamics.getDimension(); b++) {
				if (a!=b)
					out.print("migrationEvents." + mapped.dynamics.getStringStateValue(a) + "_to_"+ mapped.dynamics.getStringStateValue(b) + "\t");				
			}
		}
		out.print("nrMigrationEvents\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		mapped.calcForLogging(sample);
		// compute the number of events
		counts = new int[mapped.dynamics.getDimension()][mapped.dynamics.getDimension()];
		getCounts(mapped.getRoot());
		int totNr = 0;
		for (int a=0; a < mapped.dynamics.getDimension(); a++) {
			for (int b = 0; b < mapped.dynamics.getDimension(); b++) {
				if (a!=b) {
					out.print(counts[a][b] + "\t");
					totNr +=counts[a][b];
				}
			}
		}		
		out.print(totNr + "\t");	
		
	}


	private void getCounts(Node node) {
		if (node.getChildCount()==2) {
			getCounts(node.getLeft());
			getCounts(node.getRight());
		}else if (node.getChildCount()==1) {
			// add count
			counts[(int) node.getMetaData("location")][(int) node.getChild(0).getMetaData("location")]++;	
			getCounts(node.getChild(0));
		}		
	}

	@Override
	public void close(PrintStream out) {
	}

}
