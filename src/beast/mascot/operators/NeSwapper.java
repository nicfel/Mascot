package beast.mascot.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

public class NeSwapper extends Operator {
    public Input<List<RealParameter>> logNeInput = new Input<>(
    		"logNe", "input of the log effective population sizes", new ArrayList<>());
    public Input<RealParameter> migrationInput = new Input<>(
    		"migration", "input of the log effective population sizes");

    int length;
    int dim;
    RealParameter mig;
    
	int[][] dirs; 

    
	@Override
	public void initAndValidate() {
		dim = logNeInput.get().size();
		length = logNeInput.get().get(0).getDimension();
		mig = migrationInput.get();
		for (int i = 0; i < dim; i++) {
			if (logNeInput.get().get(0).getDimension()!=length)
				throw new IllegalArgumentException("all input paramter have to have the same dimension");
		}	
		dirs = new int[dim][dim];
		
		int c = 0;
		for (int a = 0; a < dim; a++) {
			for (int b = 0; b < dim; b++) {
				if (a!=b) {
					dirs[a][b] = c;
					c++;
				}					
			}			
		}
		
	}

	@Override
	public double proposal() {
		int nrSpots = Randomizer.nextInt(length)+1;
//		int nrSpots = 1;
		
		int startSpot = 0;
		if (nrSpots!=length) {
			startSpot= Randomizer.nextInt(length-nrSpots+1);
		}
		
		int i = Randomizer.nextInt(dim);
		int j = Randomizer.nextInt(dim);
		
		while (i == j)
			j = Randomizer.nextInt(dim);	
					
		for (int a = 0; a < nrSpots; a++) {
			int index = a+startSpot;
			double val = logNeInput.get().get(i).getArrayValue(index);
			logNeInput.get().get(i).setValue(index, logNeInput.get().get(j).getArrayValue(index));
			logNeInput.get().get(j).setValue(index, val);
		}	
				
		//
//		List<Integer> froms = new ArrayList<>();
//		List<Integer> tos = new ArrayList<>();
		
		for (int a = 0; a < dim; a++) {
//			System.out.println(Arrays.toString(dirs[a]));
			if (i!=a) {
				if (j==a) {
//					froms.add(dirs[i][a]);
//					tos.add(dirs[j][i]);
					if (Randomizer.nextBoolean()) {
						mig.swap(dirs[i][a], dirs[j][i]);
					}
				}else {
//					froms.add(dirs[i][a]);
//					tos.add(dirs[j][a]);

					if (Randomizer.nextBoolean()) {
						mig.swap(dirs[i][a], dirs[j][a]);
					}
				}
			}
		}		
		
		return 0;
	}    
    
    

}
