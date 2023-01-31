package mascot.dynamics;



import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.Dynamics;
import mascot.dynamics.RateShifts;
import mascot.parameterdynamics.NeDynamicsList;


@Description("Wrapper class that takes parametric and non parametric dynamics as input.")
public class DynamicEffectivePopulationSizesBSSVS extends Dynamics implements Loggable {	
    
    
    public Input<NeDynamicsList> parametricFunctionInput = new Input<>(
    		"NeDynamics", "input of the log effective population sizes", Validate.REQUIRED);    
    
    public Input<RealParameter> f_mInput = new Input<>(
    		"forwardsMigration", "input of backwards in time migration rates", Validate.REQUIRED);  
    
    public Input<RateShifts> rateShiftsInput = new Input<>(
    		"rateShifts", "timing of the rate shifts", Validate.REQUIRED);   
    
    public Input<BooleanParameter> indicatorInput = new Input<>(
    		"indicators", "input of backwards in time migration rates", Validate.REQUIRED);    

	public Input<Double> maxRateInput = new Input<>(
			"maxRate", "maximum rate used for integration", Double.POSITIVE_INFINITY);
    

	double[] intTimes;	
	 
	int firstlargerzero;

	boolean isForward = false;
	
	RealParameter migration;
	
	NeDynamicsList parametricFunction;
	
	int[][] dirs;
	
		
    @Override
    public void initAndValidate() {
    	
    	super.initAndValidate();
    	
    	if (dimensionInput.get()<1)
    		dimensionInput.set(getNrTypes());
    	

		migration = f_mInput.get();
		isForward = true;

    	// check the dimension of the migration rates
    	if (dimensionInput.get()*(dimensionInput.get()-1)!=migration.getDimension() || fromBeautiInput.get()){
    		System.err.println("the dimension of " + migration.getID() + " is set to " + (dimensionInput.get()*(dimensionInput.get()-1)));
    		migration.setDimension(dimensionInput.get()*(dimensionInput.get()-1));
    	}
    	
    	parametricFunction = parametricFunctionInput.get();
    	    	
    	// check that all the inputs have the correct dimension
    	if (dimensionInput.get()!=parametricFunction.size()){
    		throw new IllegalArgumentException("the number of parametric functions given as Input does not match the dimension/number of states");
    	}
    	
		//get the first non zero element
		firstlargerzero = 0;
		for (int i=0; i < rateShiftsInput.get().getDimension(); i++) {
			if (rateShiftsInput.get().getValue(i)>0) {
				firstlargerzero=i;
				break;				
			}
		}
		// initialize the intervals
		intTimes = new double[rateShiftsInput.get().getDimension()-firstlargerzero];
		for (int i=0; i < intTimes.length; i++) {
			if (i==0) {
				intTimes[i] = rateShiftsInput.get().getValue(i+firstlargerzero);						
			}
			else {
				intTimes[i] = rateShiftsInput.get().getValue(i+firstlargerzero)-rateShiftsInput.get().getValue(i-1+firstlargerzero);
			}
		}
		
		// set the Ne dimension for all skygrid dynamics.
		for (int i = 0; i < parametricFunction.size(); i++)
			parametricFunction.get(i).setNrIntervals(intTimes.length);
		
    	// Set the dimension of the indicator variables
		if (indicatorInput.get().getDimension() != migration.getDimension())
			indicatorInput.get().setDimension(migration.getDimension());
		
    	hasIndicators = true;
    	
		dirs = new int[dimensionInput.get()][dimensionInput.get()];
		
		int c = 0;
		for (int a = 0; a < dimensionInput.get(); a++) {
			for (int b = 0; b < dimensionInput.get(); b++) {
				if (a!=b) {
					dirs[a][b] = c;
					c++;
				}					
			}			
		}

    }

    /**
     * Returns the time to the next interval.
     */
    @Override
    public double getInterval(int i) {
    	if (i >= intTimes.length){
     		return Double.POSITIVE_INFINITY;
     	}else{
			return intTimes[i];
     	}
    }   

    @Override
    public double[] getIntervals() {
    	return intTimes;
    }
    
	@Override
    public double[] getCoalescentRate(int i){
		int intervalNr;
    	if (i >= rateShiftsInput.get().getDimension()-firstlargerzero-1)
    		intervalNr = rateShiftsInput.get().getDimension()-2;
    	else
    		intervalNr = i + firstlargerzero;

    	// get the current time as the midpoint of the current intervals
    	double t = rateShiftsInput.get().getIntervalMidpoint(intervalNr);   

		double[] coal = new double[getDimension()];
		for (int j = 0; j < coal.length; j++){
			if (parametricFunction.get(j).isTime)
				coal[j] = Math.min(maxRateInput.get(), 1/parametricFunction.get(j).getNeTime(t));
			else
				coal[j] = Math.min(maxRateInput.get(), 1/parametricFunction.get(j).getNeInterval(intervalNr));
		}
		return coal;
    }
    
	@Override    
    public double[] getBackwardsMigration(int i){
		int intervalNr;
    	if (i >= rateShiftsInput.get().getDimension()-firstlargerzero-1)
    		intervalNr = rateShiftsInput.get().getDimension()-2;
    	else
    		intervalNr = i + firstlargerzero;
		
    	int dim = dimensionInput.get();
    	int start = intervalNr*dim*(dim-1);
    	
		double[] m = new double[dim * dim];
		
		double[] NeInt = new double[dim];
    	double t = rateShiftsInput.get().getIntervalMidpoint(intervalNr);   	
		
		for (int j = 0; j < NeInt.length; j++){
			if (parametricFunction.get(j).isTime)
				NeInt[j] = parametricFunction.get(j).getNeTime(t);
			else
				NeInt[j] = parametricFunction.get(j).getNeInterval(intervalNr);
		}

		
		for (int a = 0; a < dim; a++){
			for (int b = 0; b < dim; b++){
				if (a!=b){
    				if (indicatorInput.get().getArrayValue(dirs[b][a])>0.5)
						m[a * dim + b] = Math.min(maxRateInput.get(), NeInt[b]*migration.getArrayValue(dirs[b][a])/NeInt[a]);
    				else
						m[a * dim + b] = 0.0;
				}
			}
		}	
		return m;
    }

//	@Override
//	public void recalculate() {
//		// TODO Auto-generated method stub
//		
//	}    
//	
//	@Override
//	public void init(PrintStream out) {
//		for (int j = 0; j < dimensionInput.get(); j++){
//			for (int i = 0; i < rateShifts.getDimension(); i++){
//				out.print(String.format("Ne.%d.%d\t", j,i));
//			}			
//		}
//	}
//
//	@Override
//	public void log(long sample, PrintStream out) {
//		for (int j = 0; j < dimensionInput.get(); j++){
//			for (int i = 0; i < rateShifts.getDimension(); i++){
//				out.print(Math.exp(Ne.getArrayValue(dimensionInput.get()*i+j))*NeMean.getValue() + "\t");
//			}			
//		}
//	}
	
	
	@Override
	public int[] getIndicators(int i) {
		int nrTrue = 0;
		for (int j = 0; j < indicatorInput.get().getDimension(); j++)
			if (indicatorInput.get().getArrayValue(j) > 0.5)
				nrTrue++;		

		int[] m = new int[nrTrue*2];
    	int c = 0;
    	int mi = 0;
    	
    	
		if (isForward){
			
	    	for (int a = 0; a < dimensionInput.get(); a++){
	    		for (int b = 0; b < dimensionInput.get(); b++){
	    			if (a!=b){
	    				if (indicatorInput.get().getArrayValue(dirs[b][a])>0.5){
	    					m[mi * 2 + 0] = a;
	    					m[mi * 2 + 1] = b;
	    					mi++;
	    				}
	    				c++;
	    			}
	    		}
	    	} 


			
		}else{
//			int c1=0,c2=0;
//			start = 0;
//			for (int a = 0; a < dim; a++){
//				for (int b = 0; b < dim; b++){
//					if (a!=b){
//	    				if (indicatorInput.get().getArrayValue(c2)>0.5)
//							m[c1] = migration.getArrayValue(start+c2);
//	    				else
//							m[c1] = 0;
//						c2++;
//					}
//					c1++;
//				}
//			}
		}	
		return m;

    	
		
		    	
//    	if (isBackwardsMigration){
//    		if (migrationType == MigrationType.asymmetric){
//		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
//		    		for (int b = 0; b < NeInput.get().getDimension(); b++){
//		    			if (a!=b){
//		    				if (indicatorInput.get().getArrayValue(c)>0.5){
//		    					m[mi * 2 + 0] = a;
//		    					m[mi * 2 + 1] = b;
//		    					mi++;
//		    				}
//		    				c++;
//		    			}
//		    		}
//		    	}
//    		}else{
//		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
//		    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
//		    			if (a!=b){
//		    				if (indicatorInput.get().getArrayValue(c)>0.5){
//		    					m[mi * 2 + 0] = a;
//		    					m[mi * 2 + 1] = b;
//		    					mi++;
//		    					m[mi * 2 + 0] = b;
//		    					m[mi * 2 + 1] = a;
//		    					mi++;
//		    				}
//		    				c++;
//		    			}
//		    		}
//		    	}    			
//    		}
//    	}else{
//    		if (migrationType == MigrationType.asymmetric){
//		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
//		    		for (int b = 0; b < NeInput.get().getDimension(); b++){
//		    			if (a!=b){
//		    				if (indicatorInput.get().getArrayValue(c)>0.5){
//		    					m[mi * 2 + 0] = a;
//		    					m[mi * 2 + 1] = b;
//		    					mi++;
//		    				}
//		    				c++;
//		    			}
//		    		}
//		    	} 
//    		}else{
//		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
//		    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
//		    			if (a!=b){
//		    				if (indicatorInput.get().getArrayValue(c)>0.5){
//		    					m[mi * 2 + 0] = a;
//		    					m[mi * 2 + 1] = b;
//		    					mi++;
//		    					m[mi * 2 + 0] = b;
//		    					m[mi * 2 + 1] = a;
//		    					mi++;
//		    				}
//		    				c++;
//		    			}
//		    		}
//		    	}    			
//    		}
//    	}       	
//    	
//    	return m;  	
	}

	

    @Override
    public int getEpochCount() {
    	return rateShiftsInput.get().getCorrectedDimension();
    }

	@Override
	public boolean intervalIsDirty(int j) {
		
		for (int i = 0; i < parametricFunction.size(); i++)
			if(parametricFunction.get(i).isDirty())
				return true;
				
		for (int i = 0; i < migration.getDimension(); i++)
			if(migration.isDirty(i))
				return true;

		if(rateShiftsInput.get().somethingIsDirty())
			return true;

		return false;
	}

	@Override
	public void init(PrintStream out) {
		if (isForward){
			for (int a = 0; a < dimensionInput.get(); a++) {
				for (int b = 0; b < intTimes.length; b++) {
					out.print(String.format("Ne_%s.%d\t", getStringStateValue(a), b));
				}
			}			

			
			for (int a = 0; a < dimensionInput.get(); a++)
				for (int b = 0; b < dimensionInput.get(); b++)
					if (a!=b)
						out.print(String.format("f_%s.%s_to_%s\t", f_mInput.get().getID(), getStringStateValue(a), getStringStateValue(b)));
			}			
	}
	
	@Override
	public void log(long sample, PrintStream out) {
		int c=0;
		if (isForward){
			for (int a = 0; a < dimensionInput.get(); a++) {
				for (int b = 0; b < intTimes.length; b++) {
					double t = rateShiftsInput.get().getIntervalMidpoint(b);
					double ne;
					if (parametricFunction.get(a).isTime)
						ne = parametricFunction.get(a).getNeTime(t);
					else
						ne = parametricFunction.get(a).getNeInterval(b);

					out.print(Math.log(ne) + "\t");
				}
			}
			
			
			for (int a = 0; a < dimensionInput.get(); a++) {
				for (int b = 0; b < dimensionInput.get(); b++) { 
					if (a!=b) {
						if (indicatorInput.get().getArrayValue(dirs[a][b])>0.5)
							out.print(f_mInput.get().getArrayValue(dirs[a][b]) + "\t");
						else
							out.print(0.0 + "\t");
						
						c++;
					}
				}	
			}
		}			
		
	}


	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void recalculate() {
		// TODO Auto-generated method stub
		
	}
	
}