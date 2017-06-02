package beast.mascot.ode;

import java.util.Arrays;

public class RungeKutta2order extends MascotSeperatedDifferentialEquation {

	public RungeKutta2order(double[][] migration_rates, double[] coalescent_rates, int lineages, int states) {
		super(migration_rates, coalescent_rates, lineages, states);
	}
	
	public void calculateValues(double duration, double[] p){
//		System.out.println(" ");
		double remainingTime = duration;
		double[] k1 = new double[p.length];
		double[] k2 = new double[p.length];
		double[] p_new = new double[p.length];
		while (remainingTime > 0){ 
			double scaleTs = Integer.MAX_VALUE;
			while(scaleTs!=0.0){
				computeDerivatives(p, k1);			
				// calculate an initial condition for the time step
				double timeStep = getTimeStep(remainingTime, p,  k1);
				remainingTime -= timeStep;			
				multiplyAdd(p_new, p, k1, timeStep);
				computeDerivatives(p_new, k2);
			}
			
			updateP(p, k1, k2, timeStep);			
		}
	}
	
	private double newTS (double[] k1, double[] k2){
		double diff = 0;
		for (int i = 0; i < k1.length-1; i++){
			double new_diff = Math.abs(k1[i] -  k2[i]);
			if (new_diff>diff)
				diff = new_diff;		
		}
		
		if (diff>0.1) return 0.1/diff;
		else return 0.0;
			
	}
	
	private double getTimeStep (double timeStep, double[] p, double[] pDot){
		double max_ratio = 1.0;
		for (int i = 0; i < (p.length-1); i++){
			double new_val = p[i] + pDot[i]*timeStep;
			double lower = p[i]*0.99;
			double upper = (1-p[i])*0.99 + p[i];
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
		return max_ratio*timeStep;		
	}
	
	private void multiplyAdd(double[] p_new, double[] p, double[] pDot, double ts){
		for (int i = 0; i < p.length; i++)
			p_new[i] = p[i] + pDot[i]*ts;	
	}

	
	private void updateP(double[] p, double[] k1, double[] k2, double[] k3, double[] k4, double timeStep){
		double new_val;
		for (int i = 0; i < p.length; i++){
			new_val = p[i] + timeStep*(k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6);
			if (new_val < 0 && i < (p.length-1)){
				System.out.println("dgkljdgjksldlksjkldfs");
				System.out.println(timeStep);
				System.out.println(" old: " + p[i] + " k1: " + k1[i] + " k2: " + k2[i] + " k3: " + k3[i] + " k4: " + k4[i]);
				System.out.println(new_val);
				System.exit(0);
			}else{
				p[i] = new_val;
			}
		}
	}
	
}


