package beast.mascot.ode;

import java.util.Arrays;

public class RungeKutta extends MascotSeperatedDifferentialEquation {

	public RungeKutta(double[][] migration_rates, double[] coalescent_rates, int lineages, int states) {
		super(migration_rates, coalescent_rates, lineages, states);
	}
	
	public void calculateValues(double duration, double[] p){
//		System.out.println(" ");
		double remainingTime = duration;
		double[] k1 = new double[p.length];
		double[] k2 = new double[p.length];
		double[] k3 = new double[p.length];
		double[] k4 = new double[p.length];
		double[] p_new = new double[p.length];
		while (remainingTime > 0){ 
			computeDerivatives(p, k1);
//			System.out.println("p: " + Arrays.toString(p));
//			System.out.println("k1: " + Arrays.toString(k1));
			// calculate a good time step based on the euler calculation
			double timeStep = getTimeStep(remainingTime, p,  k1);
			timeStep /= 100;
			remainingTime -= timeStep;			
			multiplyAdd(p_new, p, k1, timeStep/2);
//			System.out.println("p1: " + Arrays.toString(p_new));
			computeDerivatives(p_new, k2);
			multiplyAdd(p_new, p, k2, timeStep/2);
//			System.out.println("p2: " + Arrays.toString(p_new));
			computeDerivatives(p_new, k3);
			multiplyAdd(p_new, p, k3, timeStep);
//			System.out.println("p3: " + Arrays.toString(p_new));
			computeDerivatives(p_new, k4);
			updateP(p, k1, k2, k3, k4, timeStep);			
//			System.out.println(timeStep);
//			System.out.println("k  " + Arrays.toString(k1));
		}
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


