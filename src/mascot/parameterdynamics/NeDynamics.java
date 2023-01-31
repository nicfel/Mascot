package mascot.parameterdynamics;

import beast.base.inference.CalculationNode;

public abstract class NeDynamics extends CalculationNode {
	
	public boolean isTime;
	
	@Override
	public void initAndValidate() {
	}

	
    /**
     * recalculate the dynamics
     */
    public void recalculate() {};

    /**
     * returns the effective population size at time t
     * @param t
     * @return
     */
    public double getNeTime(double t) {
		throw new IllegalArgumentException("Function not implemented. Class of parametric function not correctly recognized");
    }
    

    /**
     * returns the effective population size at interval i
     * @param i
     * @return
     */
    public double getNeInterval(int i) {
		throw new IllegalArgumentException("Function not implemented. Class of parametric function not correctly recognized");
    }
    
    public void setNrIntervals(int intervals) {}

    
	public boolean isDirty() {
		return true;
	};
    
    
}
