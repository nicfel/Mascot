package beast.mascot.dynamics;


import java.io.PrintStream;
import java.util.Arrays;

import org.jblas.DoubleMatrix;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;


@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class Constant extends Dynamics implements Loggable  {

    public Input<RealParameter> NeInput = new Input<>("Ne", "input of effective population sizes", Validate.REQUIRED);    
    public Input<RealParameter> b_mInput = new Input<>("backwardsMigration", "input of backwards in time migration rates");    
    public Input<RealParameter> f_mInput = new Input<>("forwardsMigration", "input of backwards in time migration rates", Validate.XOR, b_mInput);    

	private boolean isBackwardsMigration;
	
	private enum MigrationType {
	    symmetric, asymmetric 
	}
	
	MigrationType migrationType;

    

    @Override
    public void initAndValidate() {
    	super.initAndValidate();
    	
    	
    	if (b_mInput.get()!=null)
    		isBackwardsMigration = true;
    	else
    		isBackwardsMigration = false;
    	
    	int migDim = NeInput.get().getDimension()*(NeInput.get().getDimension()-1);
    	
    	if (isBackwardsMigration){
    		if (migDim == b_mInput.get().getDimension())
    			migrationType = MigrationType.asymmetric;
    		else if ((int) migDim/2 == b_mInput.get().getDimension())
    			migrationType = MigrationType.symmetric;
    		else
    			System.err.println("wrong number of migraiton or Ne elements");
    	}
    	if (!isBackwardsMigration){
    		if (migDim == f_mInput.get().getDimension()){
    			migrationType = MigrationType.asymmetric;
    		}else if ((int) migDim/2 == f_mInput.get().getDimension()){
    			migrationType = MigrationType.symmetric;
    		}else{
    			System.err.println("wrong number of migration or Ne elements");
    			System.err.println("check the dimension of the migration and Ne parameter in the xml");
//    			System.exit(0);
    		}
    	}
   	
    }

    /**
     * Returns the time to the next interval.
     */
    public double getInterval(int i) {
    	return Double.POSITIVE_INFINITY;
    }   
    
    public boolean intervalIsDirty(int i){
    	boolean intervalIsDirty = false;  	    	
    	
    	for (int j = 0; j < dimensionInput.get(); j++)
    		if (NeInput.get().isDirty(j))
    			intervalIsDirty = true;  
    
    	if (isBackwardsMigration){
			for (int j = 0; j < b_mInput.get().getDimension(); j++)
				if (b_mInput.get().isDirty(j))
					intervalIsDirty = true; 
    	}else{
			for (int j = 0; j < f_mInput.get().getDimension(); j++)
				if (f_mInput.get().isDirty(j))
					intervalIsDirty = true;    		
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
    public double[][] getBackwardsMigration(int i){
    	double[][] m = new double[NeInput.get().getDimension()][NeInput.get().getDimension()];
    	
    	if (isBackwardsMigration){
    		if (migrationType == MigrationType.asymmetric){
		    	int c = 0;
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = 0; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				m[a][b] = b_mInput.get().getArrayValue(c);
		    				c++;
		    			}
		    		}
		    	}
    		}else{
		    	int c = 0;
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				m[a][b] = b_mInput.get().getArrayValue(c);
		    				m[b][a] = b_mInput.get().getArrayValue(c);
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
		    				m[a][b] = f_mInput.get().getArrayValue(c)
		    						*NeInput.get().getArrayValue(b)
		    							/NeInput.get().getArrayValue(a);
		    				c++;
		    			}
		    		}
		    	} 
    		}else{
		    	int c = 0;
		    	for (int a = 0; a < NeInput.get().getDimension(); a++){
		    		for (int b = a+1; b < NeInput.get().getDimension(); b++){
		    			if (a!=b){
		    				m[a][b] = f_mInput.get().getArrayValue(c)
		    						*NeInput.get().getArrayValue(b)
		    							/NeInput.get().getArrayValue(a);
		    				m[b][a] = f_mInput.get().getArrayValue(c)
		    						*NeInput.get().getArrayValue(a)
		    							/NeInput.get().getArrayValue(b);
		    				c++;
		    			}
		    		}
		    	}    			
    		}
    	}
    	
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
	public void log(int sample, PrintStream out) {
		for (int i = 0 ; i < NeInput.get().getDimension(); i++){
			out.print(String.format("%f\t", NeInput.get().getArrayValue(i)));
		}
		
		if (migrationType == MigrationType.asymmetric){
	    	int c = 0;
	    	for (int a = 0; a < NeInput.get().getDimension(); a++){
	    		for (int b = 0; b < NeInput.get().getDimension(); b++){
	    			if (a!=b){
	    				if (isBackwardsMigration)
	    					out.print(String.format("%f\t", b_mInput.get().getArrayValue(c)));
	    				else
	    					out.print(String.format("%f\t", f_mInput.get().getArrayValue(c)));
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
	    					out.print(String.format("%f\t", b_mInput.get().getArrayValue(c)));
	    				else
	    					out.print(String.format("%f\t", f_mInput.get().getArrayValue(c)));
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



    
}