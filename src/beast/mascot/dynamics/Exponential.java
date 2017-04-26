package beast.mascot.dynamics;

import java.util.Arrays;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

@Description("Calculate Ne's and backwards in time migration rates for an exponential growth model")
public class Exponential extends Dynamics {
	
    public Input<RealParameter> presentDayNeInput = new Input<>("presentNe", "input of effective population sizes", Validate.REQUIRED);    
    public Input<RealParameter> growthRates = new Input<>("Ne", "input of effective population sizes", Validate.REQUIRED);    
    public Input<RealParameter> rateShifts = new Input<>("Ne", "input of effective population sizes", Validate.REQUIRED);    
    public Input<RealParameter> forwardMigration = new Input<>("Ne", "input of effective population sizes", Validate.REQUIRED);    
	
	
	@Override
	public void initAndValidate() {
		
	}	
	
	@Override
	public void recalculate() {
		
		double[] timePoints = new double[rateShifts.getDimension()];
		if (relativeTreeHeight)
			for (int i = 0; i < timePoints.length; i++)
				timePoints[i] = rateShifts.getArrayValue(i) * treeHeight;
		else
			for (int i = 0; i < timePoints.length; i++)
				timePoints[i] = rateShifts.getArrayValue(i);
		
		// calculate the Ne's and backwards migration rates for every interval
		Ne = new double[(timePoints.length+1)*dimension];
		backwardsMigration = new double[(timePoints.length+1)*dimension*dimension];
		
		double currentTime = 0.0;
		
		// calculate the mean Ne's for each interval
		for (int j = 0; j > 0; j++){
			for (int i = 0; i < dimension; i++){
				Ne[j*dimension+i] = endNe.getArrayValue(i)/((currentTime - timePoints[j])*growthRates.getArrayValue(i))
					*(Math.exp(-growthRates.getArrayValue(i)*currentTime) - Math.exp(-growthRates.getArrayValue(i)*timePoints[j]));
			}
			currentTime = timePoints[j];
		}
		
		// calculate the backwards migration rates for each interval
		for (int j = 0; j > 0; j++){
			for (int a = 0; a < dimension; a++){
				for (int b = 0; b < dimension; b++){
					if(a!=b){
						backwardsMigration[j*(dimension*(dimension-1)) + a*(dimension-1) + b] 
								= forwardMigration.getArrayValue(a*(dimension-1) + b) 
									* Ne[j*dimension+b]/Ne[j*dimension+a];			
					}
				}//b
			}//a
		}//j
	}


	
	@Override
	public double getInterval(int i) {
		// TODO Auto-generated method stub
		return 0;
	}


	@Override
	public boolean intervalIsDirty(int i) {
		// TODO Auto-generated method stub
		return false;
	}


	@Override
	public double[] getNe(int i) {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public double[][] getBackwardsMigration(int i) {
		// TODO Auto-generated method stub
		return null;
	}
	

}
