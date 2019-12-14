package beast.mascot.dynamics;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.util.FastMath;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.mascot.glmmodel.GlmModel;


@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
@Citation(	"Nicola F. Müller, Gytis Dudas, Tanja Stadler (2018)\n"+
			"  Inferring time-dependent migration and coalescence patterns\n" +
			"  from genetic sequence and predictor data in structured populations\n"+
			"  bioRxiv, doi: bty406, 10.1101/342329")
public class GLM extends Dynamics implements Loggable {	
    
	public Input<GlmModel> migrationGLMInput = new Input<>(
			"migrationGLM", "input of migration GLM model", Validate.REQUIRED);
    
	public Input<GlmModel> NeGLMInput = new Input<>(
			"NeGLM", "input of migration GLM model", Validate.REQUIRED);
    
    public Input<RateShifts> rateShiftsInput = new Input<>(
    		"rateShifts", "input of timings of rate shifts relative to the most recent sample", Validate.OPTIONAL);    
     
	public Input<Double> maxRateInput = new Input<>(
			"maxRate", "maximum rate used for integration", Double.POSITIVE_INFINITY);
	
	double[] intTimes;	
 
	int firstlargerzero;
	
	public GLM(){
    	
	}
	
    @Override
    public void initAndValidate() {
    	super.initAndValidate(); 
		rateShiftsInput.get().initAndValidate();
		migrationGLMInput.get().covariateListInput.get().initAndValidate();
		NeGLMInput.get().covariateListInput.get().initAndValidate();

		// set the number of intervals for the GLM models
		if (fromBeautiInput.get()) {
			// the type to trait map is needed to read in predictors
			migrationGLMInput.get().covariateListInput.get().traitToType = new HashMap<>(traitToType);
			NeGLMInput.get().covariateListInput.get().traitToType = new HashMap<>(traitToType);

			migrationGLMInput.get().covariateListInput.get().nrIntervals = rateShiftsInput.get().getDimension();
			NeGLMInput.get().covariateListInput.get().nrIntervals = rateShiftsInput.get().getDimension();
			
			migrationGLMInput.get().setNrDummy();
			NeGLMInput.get().setNrDummy();
		}
		if (rateShiftsInput.get().getDimension()>0) {	
			if (migrationGLMInput.get().covariateListInput.get().size()>0)
				migrationGLMInput.get().setNrIntervals(rateShiftsInput.get().getDimension(), dimensionInput.get(), true);
			if (NeGLMInput.get().covariateListInput.get().size()>0)
				NeGLMInput.get().setNrIntervals(rateShiftsInput.get().getDimension(), dimensionInput.get(), false);
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
    
    public boolean intervalIsDirty(int i){
		if(NeGLMInput.get().isDirty())
			return true;
		if(migrationGLMInput.get().isDirty())
			return true;
    	return false;
    }  
    

    
	@Override
    public double[] getCoalescentRate(int i){
		int intervalNr;
    	if (i >= rateShiftsInput.get().getDimension()-firstlargerzero-1)
    		intervalNr = rateShiftsInput.get().getDimension()-2;
    	else
    		intervalNr = i + firstlargerzero;
    	
    	double[] Ne = NeGLMInput.get().getRates(intervalNr);
		double[] coal = new double[Ne.length];
		for (int j = 0; j < Ne.length; j++){
			coal[j] = FastMath.min(1/Ne[j], maxRateInput.get());
		}
		return coal;
    }
    
	@Override    
    public double[] getBackwardsMigration(int i){
		int intervalNr;
    	if (i >= rateShiftsInput.get().getDimension()-firstlargerzero)
    		intervalNr = rateShiftsInput.get().getDimension()-1;
    	else
    		intervalNr = i + firstlargerzero;

    	int n = dimensionInput.get();
    	double[] m = new double[n * n];
		double[] mig = migrationGLMInput.get().getRates(intervalNr);
		double[] Ne = NeGLMInput.get().getRates(intervalNr);
		
		int c = 0;
		for (int a = 0; a < dimensionInput.get(); a++){
			for (int b = 0; b < dimensionInput.get(); b++){
				if (a!=b){
					m[b * n + a] = FastMath.min( 
							Ne[a]*mig[c]/Ne[b],
							maxRateInput.get());
					c++;
				}
			}
		}
		return m;
    }

	@Override
	public void recalculate() {
		// TODO Auto-generated method stub
		
	}    
	
	public Double[] getAllCoalescentRate() {
		Double[] coal = new Double[NeGLMInput.get().nrIntervals*NeGLMInput.get().verticalEntries];
		
		for (int i = 0; i < intTimes.length; i++){
	    	double[] Ne = NeGLMInput.get().getRates(i);
	    	for (int j = 0; j < Ne.length; j++)
	    		coal[i*NeGLMInput.get().verticalEntries + j] = 1/Ne[j];
		}
		return coal;
	}

	public Double[] getAllBackwardsMigration() {
		Double[] mig = new Double[migrationGLMInput.get().nrIntervals*migrationGLMInput.get().verticalEntries];
		
		for (int i = 0; i < intTimes.length; i++){
	    	double[] m = migrationGLMInput.get().getRates(i);
	    	for (int j = 0; j < m.length; j++)
	    		mig[i*migrationGLMInput.get().verticalEntries + j] = m[j];
		}
		return mig;
	}

	@Override
	public void init(PrintStream out) {
		for (int j = 0; j < dimensionInput.get(); j++){
			for (int i = 0; i < intTimes.length; i++){
				out.print(String.format("Ne.%d.%d\t", j,i));
			}			
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int j = 0; j < dimensionInput.get(); j++){
			for (int i = 0; i < intTimes.length; i++){
		    	double[] Ne = NeGLMInput.get().getRates(i);
				out.print(Ne[j] + "\t");
			}			
		}
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

//    @Override
//	protected boolean requiresRecalculation(){
//    	
//    	return intervalIsDirty(0);
//    }


    @Override
    public int getEpochCount() {
    	return rateShiftsInput.get().getDimension();
    }
}