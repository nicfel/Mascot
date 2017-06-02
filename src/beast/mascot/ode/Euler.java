package beast.mascot.ode;

import java.util.Arrays;

public class Euler extends MascotSeperatedDifferentialEquation {

	double epsilon = 0.01;
	
	public Euler(double[][] migration_rates, double[] coalescent_rates, int lineages, int states) {
		super(migration_rates, coalescent_rates, lineages, states);
	}
	
	public void calculateValues(double duration, double[] p){
		double[] pDot = new double[p.length];
		double[] pDotDot = new double[p.length];
		while (duration > 0){
			computeDerivatives(p, pDot);
			computeSecondDerivate(p, pDot, pDotDot);
			duration = updateP(duration, p,  pDot, pDotDot);		
		}
	}
	
	private double updateP (double duration, double[] p, double[] pDot, double[] pDotDot){
		double max_dotdot = 0.0;
		for (int i = 0; i < (p.length); i++){
			if (pDotDot[i] > max_dotdot)
				max_dotdot = pDotDot[i];
		}
		double max_step = Math.min(Math.sqrt(epsilon*2/max_dotdot), duration);

		double max_ratio = 1.0;
		for (int i = 0; i < (p.length-1); i++){
			double new_val = p[i] + pDot[i]*max_step;
			double lower = p[i]/2;
			double upper = (1-p[i])/2 + p[i];
			if (new_val < lower){
				double ratio = 	(lower)/(p[i]-new_val);
				if (ratio < max_ratio)
					max_ratio = ratio;				
			}
			if(new_val > upper){
				double ratio = (upper)/(new_val-p[i]);
				if (ratio < max_ratio)
					max_ratio = ratio;
			}
		}
		double stepSize = max_step*max_ratio;
		
		doUpdating(stepSize, p, pDot);
		duration -= stepSize;
		return duration;
		
		
	}
	
	private void doUpdating(double stepSize, double[] p, double[] pDot){
		for (int i = 0; i < (p.length); i++)
			p[i] += pDot[i]*stepSize;		
	}
	
}


