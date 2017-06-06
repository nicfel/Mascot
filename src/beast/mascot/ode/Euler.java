package beast.mascot.ode;

import java.util.Arrays;

import org.apache.commons.math4.util.FastMath;

public class Euler extends MascotSeperatedDifferentialEquation {

	double epsilon = 0.1;
	

	
	public Euler(double[][] migration_rates, double[] coalescent_rates, int lineages, int states) {
		super(migration_rates, coalescent_rates, lineages, states);
	}
	
	public void calculateValues(double duration, double[] p){
		double[] pDot = new double[p.length];
		double[] pDotDot = new double[p.length];
		while (duration > 0){
			computeDerivatives(p, pDot);
			approximateSecondDerivate(p, pDot, pDotDot);
			duration = updateP(duration, p,  pDot, pDotDot);		
		}
	}
	
	private double updateP (double duration, double[] p, double[] pDot, double[] pDotDot){
		double max_dotdot = 0.0;
		for (int i = 0; i < (p.length-1); i++){
			double tmp = FastMath.abs(pDotDot[i]);
			if (tmp > max_dotdot)
				max_dotdot = tmp;
		}
		double timeStep = FastMath.min(FastMath.sqrt(epsilon*2/max_dotdot), duration);

		for (int i = 0; i < (p.length-1); i++){
			double new_val = p[i] + pDot[i]*timeStep;
			double lower = p[i]*0.01;
			double upper = (1-p[i])*0.99 + p[i];
			
			while (new_val > upper || new_val < lower){
				timeStep *= 0.9;
				new_val = p[i] + pDot[i]*timeStep ;
			}			
		}
		doUpdating(timeStep, p, pDot);
		duration -= timeStep;
		return duration;
		
		
	}
	
	private void doUpdating(double stepSize, double[] p, double[] pDot){
		for (int i = 0; i < (p.length); i++)
			p[i] += pDot[i]*stepSize;		
	}
	
}


