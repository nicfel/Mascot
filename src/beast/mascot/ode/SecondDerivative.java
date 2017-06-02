package beast.mascot.ode;

import java.util.Arrays;

public class SecondDerivative extends MascotSeperatedDifferentialEquation {

	public SecondDerivative(double[][] migration_rates, double[] coalescent_rates, int lineages, int states) {
		super(migration_rates, coalescent_rates, lineages, states);
	}
	
	public void calculateValues(double duration, double[] p, double maxStep, double max_bound){
		double[] pDot = new double[p.length];
		double[] pDotDot = new double[p.length];
//		System.out.println(Arrays.toString(p));
//		System.out.println("......");
		while (duration > 0){
//			System.out.println(" ");
			computeDerivatives(p, pDot);
			computeSecondDerivate(p, pDot, pDotDot);
			double stepTime = updateP(Math.min(duration, maxStep), p,  pDot, pDotDot, max_bound);
			duration -= stepTime;
//			System.out.println(p[0]+p[1]+p[2]+p[3]+p[4]+p[5]);
//			System.out.println(Arrays.toString(p));
//			System.out.println(Arrays.toString(pDot));
//			System.out.println(Arrays.toString(pDotDot));
//			System.out.println(duration);
		}
//		System.exit(0);
	}
	
	private double updateP (double timeStep, double[] p, double[] pDot, double[] pDotDot, double max_bound){
		double stepSize = timeStep;
		double max_ratio = 1.0;
		for (int i = 0; i < (p.length-1); i++){
			double new_val = p[i] + pDot[i]*stepSize + pDotDot[i]*stepSize*stepSize/2;
			double lower = p[i]*(max_bound);
			double upper = (1-p[i])*(1-max_bound) + p[i];
			while (new_val > upper || new_val < lower){
				stepSize *= 0.8;
				new_val = p[i] + pDot[i]*stepSize + pDotDot[i]*stepSize*stepSize/2;
			}
		}		
		
		double new_val = p[p.length-1] + pDot[p.length-1]*stepSize + pDotDot[p.length-1]*stepSize*stepSize/2;
		
		while (new_val > 0){
			stepSize *= 0.8;
			new_val = p[p.length-1] + pDot[p.length-1]*stepSize + pDotDot[p.length-1]*stepSize*stepSize/2;
		}
	
		doUpdating(stepSize, p, pDot, pDotDot);
		
		return stepSize;
		
		
	}
	
	private void doUpdating(double stepSize, double[] p, double[] pDot, double[] pDotDot){
		for (int i = 0; i < (p.length); i++)
			p[i] += pDot[i]*stepSize + pDotDot[i]*stepSize*stepSize/2;		
	}
	
}


