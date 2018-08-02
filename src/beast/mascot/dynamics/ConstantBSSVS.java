package beast.mascot.dynamics;


import java.io.PrintStream;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;


@Description("Allows to use indicators on the migration rates " +
			 "to sample significant routes of migration")
public class ConstantBSSVS extends Dynamics implements Loggable  {

    public Input<RealParameter> NeInput = new Input<>("Ne", "input of effective population sizes", Validate.REQUIRED);    
    public Input<RealParameter> b_mInput = new Input<>("backwardsMigration", "input of backwards in time migration rates");    
    public Input<RealParameter> f_mInput = new Input<>("forwardsMigration", "input of backwards in time migration rates", Validate.XOR, b_mInput);    
    public Input<BooleanParameter> indicatorInput = new Input<>("indicators", "input of backwards in time migration rates", Validate.REQUIRED);    

    public Input<RealParameter> migrationClockInput = new Input<>("migrationClock", "input of backwards in time migration rates", Validate.REQUIRED);    

    
	private boolean isBackwardsMigration;
	
	private enum MigrationType {
	    symmetric, asymmetric 
	}
	
	MigrationType migrationType;

    

    @Override
    public void initAndValidate() {    	
    	super.initAndValidate();
    	
    	if (dimensionInput.get()<1)
    		dimensionInput.set(getNrTypes());

    	
    	if (b_mInput.get()!=null)
    		isBackwardsMigration = true;
    	else
    		isBackwardsMigration = false;
    	
    	
    	if (dimensionInput.get()!=NeInput.get().getDimension()){
    		System.err.println("the dimension of " + NeInput.get().getID() + " is set to " + dimensionInput.get());
    		NeInput.get().setDimension(dimensionInput.get());
    	}

    	
    	int migDim = NeInput.get().getDimension()*(NeInput.get().getDimension()-1);
    	
    	if (isBackwardsMigration){
    		if (migDim == b_mInput.get().getDimension()){
    			migrationType = MigrationType.asymmetric;
    		}else if ((int) migDim/2 == b_mInput.get().getDimension()){
    			migrationType = MigrationType.symmetric;
    		}else{
    			migrationType = MigrationType.asymmetric;
    			System.err.println("Wrong number of migration elements, assume asymmetric migration:");
        		System.err.println("the dimension of " + b_mInput.get().getID() + " is set to " + migDim);
    			b_mInput.get().setDimension(migDim);       		
    		}
    	}
    	if (!isBackwardsMigration){
    		if (migDim == f_mInput.get().getDimension()){
    			migrationType = MigrationType.asymmetric;
    		}else if ((int) migDim/2 == f_mInput.get().getDimension()){
    			migrationType = MigrationType.symmetric;
    		}else{
    			migrationType = MigrationType.asymmetric;
    			System.err.println("Wrong number of migration elements, assume asymmetric migration:");
        		System.err.println("the dimension of " + f_mInput.get().getID() + " is set to " + migDim);
        		f_mInput.get().setDimension(migDim);       		
    		}
    	}
    	
    	// Set the dimension of the indicator variables
    	if (isBackwardsMigration)
    		indicatorInput.get().setDimension(b_mInput.get().getDimension());
		else
			indicatorInput.get().setDimension(f_mInput.get().getDimension());
    	
    	hasIndicators = true;
   	
    }

    /**
     * Returns the time to the next interval.
     */
    @Override
    public double getInterval(int i) {
    	return Double.POSITIVE_INFINITY;
    }   

    @Override
    public double [] getIntervals() {
		return new double[]{Double.POSITIVE_INFINITY};
	}

    public boolean intervalIsDirty(int i){
    	boolean intervalIsDirty = false;  
    	
    	if (migrationClockInput.get().isDirty(0))
			intervalIsDirty = true;

    	
    	for (int j = 0; j < dimensionInput.get(); j++)
    		if (NeInput.get().isDirty(j))
    			intervalIsDirty = true;
    	
    	if (isBackwardsMigration){
			for (int j = 0; j < b_mInput.get().getDimension(); j++){
				if (indicatorInput.get().isDirty(j))
					intervalIsDirty = true; 
				if (b_mInput.get().isDirty(j))
					intervalIsDirty = true; 
			}
    	}else{
			for (int j = 0; j < f_mInput.get().getDimension(); j++){
				if (indicatorInput.get().isDirty(j))
					intervalIsDirty = true; 
				if (f_mInput.get().isDirty(j))
					intervalIsDirty = true;    		
			}
    	}		
    	return intervalIsDirty;

    }   
 
	@Override
    public double[] getCoalescentRate(int i){
    	double[] coal = new double[NeInput.get().getDimension()];
    	int c = 0;
    	for (int j = 0; j < NeInput.get().getDimension(); j++){
    		coal[c] = 1/(2*NeInput.get().getArrayValue(j));
    		c++;
    	}
    	return coal;

    }
    
	@Override    
    public double[] getBackwardsMigration(int i){
		int n = NeInput.get().getDimension();
    	double[] m = new double[n * n];
//    	System.out.println(b_mInput.get());
    	
    	if (isBackwardsMigration){
    		if (migrationType == MigrationType.asymmetric){
		    	int c = 0;
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = 0; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				if (indicatorInput.get().getArrayValue(c)>0.5)
		    					m[a * n + b] = migrationClockInput.get().getValue()*b_mInput.get().getArrayValue(c);
		    				else
		    					m[a * n + b] = 0.0;
		    				c++;
		    			}
		    		}
		    	}
    		}else{
		    	int c = 0;
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				if (indicatorInput.get().getArrayValue(c)>0.5){
			    				m[a * n + b] = migrationClockInput.get().getValue()*b_mInput.get().getArrayValue(c);
			    				m[b * n + a] = migrationClockInput.get().getValue()*b_mInput.get().getArrayValue(c);
		    				}else{
		    					m[a * n + b] = 0.0;
		    					m[b * n + a] = 0.0;
		    				}
		    				c++;
		    			}
		    		}
		    	}    			
    		}
    	}else{
    		if (migrationType == MigrationType.asymmetric){
		    	int c = 0;
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = 0; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				if (indicatorInput.get().getArrayValue(c)>0.5)
		    					m[a * n + b] = migrationClockInput.get().getValue()*f_mInput.get().getArrayValue(c)
			    						*NeInput.get().getArrayValue(b)
			    							/NeInput.get().getArrayValue(a);
		    				else
		    					m[a * n + b] = 0.0;
		    				c++;
		    			}
		    		}
		    	} 
    		}else{
		    	int c = 0;
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				if (indicatorInput.get().getArrayValue(c)>0.5){
			    				m[a * n + b] = migrationClockInput.get().getValue()*f_mInput.get().getArrayValue(c)
			    						*NeInput.get().getArrayValue(b)
			    							/NeInput.get().getArrayValue(a);
			    				m[b * n + a] = migrationClockInput.get().getValue()*f_mInput.get().getArrayValue(c)
			    						*NeInput.get().getArrayValue(a)
			    							/NeInput.get().getArrayValue(b);
		    				}else{
		    					m[a * n + b] = 0.0;
		    					m[b * n + a] = 0.0;
		    				}
		    				c++;
		    			}
		    		}
		    	}    			
    		}
    	}
//    	for (int a = 0; a < getDimension(); a++){
//    		for (int b = 0; b < getDimension()-1; b++){
//    			if (m[a][b]==0)
//        			System.out.print("0, ");
//    			else
//    				System.out.print(m[a][b] + ", ");
//    		}
//			if (m[a][getDimension()-1]==0)
//    			System.out.print("0;");
//			else
//				System.out.print(m[a][getDimension()-1] + ";");
//    		System.out.print("\n");
//    	}
//    	System.out.println();
    	return m;  	
    }    

	@Override 
	public int getDimension() {
		return NeInput.get().getDimension();
	}
	
	@Override
	public void recalculate() {
	}

	@Override
	public void init(PrintStream out) {
		for (int i = 0 ; i < NeInput.get().getDimension(); i++){
			out.print(String.format("Ne.%s\t", getStringStateValue(i)));
		}
		if (migrationType == MigrationType.asymmetric){
	    	int c = 0;
	    	for (int a = 0; a < NeInput.get().getDimension(); a++){
	    		for (int b = 0; b < NeInput.get().getDimension(); b++){
	    			if (a!=b){
	    				if (isBackwardsMigration)
	    					out.print(String.format("b_migration.%s_to_%s\t", getStringStateValue(a), getStringStateValue(b)));
	    				else
	    					out.print(String.format("f_migration.%s_to_%s\t", getStringStateValue(b), getStringStateValue(a)));	    					
	    				c++;
	    			}
	    		}
	    	} 
		}
		else{
			int c = 0;
	    	for (int a = 0; a < NeInput.get().getDimension(); a++){
	    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
	    			if (a!=b){
	    				if (isBackwardsMigration)
	    					out.print(String.format("b_migration.%s_and_%s\t", getStringStateValue(a), getStringStateValue(b)));
	    				else
	    					out.print(String.format("f_migration.%s_and_%s\t", getStringStateValue(b), getStringStateValue(a)));	    					
	    				c++;
	    			}
	    		}
	    	} 
		}

	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0 ; i < NeInput.get().getDimension(); i++){
			out.print(String.format("%f\t", NeInput.get().getArrayValue(i)));
		}
		
		if (migrationType == MigrationType.asymmetric){
	    	int c = 0;
	    	for (int a = 0; a < NeInput.get().getDimension(); a++){
	    		for (int b = 0; b < NeInput.get().getDimension(); b++){
	    			if (a!=b){
	    				if (isBackwardsMigration){
		    				if (indicatorInput.get().getArrayValue(c)>0.5)
		    					out.print(String.format("%f\t", migrationClockInput.get().getValue()*b_mInput.get().getArrayValue(c)));
		    				else
		    					out.print("0.0\t");
	    				}else{
		    				if (indicatorInput.get().getArrayValue(c)>0.5)
		    					out.print(String.format("%f\t", migrationClockInput.get().getValue()*f_mInput.get().getArrayValue(c)));
		    				else
		    					out.print("0.0\t");
	    				}
	    				c++;
	    			}
	    		}
	    	} 
		}
		else{
			int c = 0;
	    	for (int a = 0; a < NeInput.get().getDimension(); a++){
	    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
	    			if (a!=b){
	    				if (isBackwardsMigration){
		    				if (indicatorInput.get().getArrayValue(c)>0.5)
		    					out.print(String.format("%f\t", migrationClockInput.get().getValue()*b_mInput.get().getArrayValue(c)));
		    				else
		    					out.print("0.0\t");
	    				}else{
		    				if (indicatorInput.get().getArrayValue(c)>0.5)
		    					out.print(String.format("%f\t", migrationClockInput.get().getValue()*f_mInput.get().getArrayValue(c)));
		    				else
		    					out.print("0.0\t");
	    				}
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
	public int[] getIndicators(int i) {
		int nrTrue = 0;
		for (int j = 0; j < indicatorInput.get().getDimension(); j++)
			if (indicatorInput.get().getArrayValue(j) > 0.5)
				nrTrue++;
		
		if (migrationType == MigrationType.symmetric)
			nrTrue*=2;
		
		int[] m = new int[nrTrue*2];
    	int c = 0;
    	int mi = 0;
		
		    	
    	if (isBackwardsMigration){
    		if (migrationType == MigrationType.asymmetric){
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = 0; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				if (indicatorInput.get().getArrayValue(c)>0.5){
		    					m[mi * nrTrue + 0] = a;
		    					m[mi * nrTrue + 1] = b;
		    					mi++;
		    				}
		    				c++;
		    			}
		    		}
		    	}
    		}else{
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				if (indicatorInput.get().getArrayValue(c)>0.5){
		    					m[mi * nrTrue + 0] = a;
		    					m[mi * nrTrue + 1] = b;
		    					mi++;
		    					m[mi * nrTrue + 0] = b;
		    					m[mi * nrTrue + 1] = a;
		    					mi++;
		    				}
		    				c++;
		    			}
		    		}
		    	}    			
    		}
    	}else{
    		if (migrationType == MigrationType.asymmetric){
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = 0; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				if (indicatorInput.get().getArrayValue(c)>0.5){
		    					m[mi * nrTrue + 0] = a;
		    					m[mi * nrTrue + 1] = b;
		    					mi++;
		    				}
		    				c++;
		    			}
		    		}
		    	} 
    		}else{
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				if (indicatorInput.get().getArrayValue(c)>0.5){
		    					m[mi * nrTrue + 0] = a;
		    					m[mi * nrTrue + 1] = b;
		    					mi++;
		    					m[mi * nrTrue + 0] = b;
		    					m[mi * nrTrue + 1] = a;
		    					mi++;
		    				}
		    				c++;
		    			}
		    		}
		    	}    			
    		}
    	}
    	
    	return m;  	
	}

//    @Override
//	protected boolean requiresRecalculation(){
//    	return intervalIsDirty(0);
//    }
//
    @Override
    public int getEpochCount() {return 1;}

    
}