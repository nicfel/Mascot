package beast.mascot.parameterDynamics;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

@Description("dummy class for Beauti")
public class NeDynamics extends CalculationNode {
	public Input<EffectivePopulationSizeDynamics> dynamicsInput = new Input<>("dyn", 
			"input of effective populaiton size dynamics", Validate.REQUIRED);
	
	public boolean isTime;
	
	@Override
	public void initAndValidate() {
		dynamicsInput.get().initAndValidate();
	}

	
    /**
     * recalculate the dynamics
     */
    public void recalculate() {
		dynamicsInput.get().recalculate();

    };

    /**
     * returns the effective population size at time t
     * @param t
     * @return
     */
    public double getNeTime(double t) {
		return dynamicsInput.get().getNeTime(t);
    }
    

    /**
     * returns the effective population size at interval i
     * @param i
     * @return
     */
    public double getNeInterval(int i) {
		return dynamicsInput.get().getNeInterval(i);
    }
    
    public void setNrIntervals(int intervals) {
    	dynamicsInput.get().setNrIntervals(intervals);
    }

    
	public boolean isDirty() {
    	return dynamicsInput.get().isDirty();
	};
    

}
