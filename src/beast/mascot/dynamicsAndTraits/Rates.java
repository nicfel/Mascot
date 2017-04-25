package beast.mascot.dynamicsAndTraits;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math4.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math4.ode.FirstOrderIntegrator;
import org.apache.commons.math4.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math4.ode.nonstiff.HighamHall54Integrator;
import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.mascot.ode.SIR;

/**
 * @author Nicola Felix Mueller 
 */
@Description("Calculates the probability of particels under an SIR model")
public class Rates extends CalculationNode implements Loggable {
	// Input in case of constant rates
    public Input<RealParameter> constantCoalescentRatesInput = new Input<>("constantCoalescentRates", "input of the effective populationSize");
    public Input<RealParameter> constantMigrationRatesInput = new Input<>("constantMigrationRates", "input of the effective populationSize");
    public Input<BooleanParameter> constantMigrationIndicatorInput = new Input<>("constantMigrationIndicator", "indicates if a migration rate is zero or non zero");     
    
    // input in case of having an SIR Type model
    public Input<RealParameter> transmissionRatesInput = new Input<>("transmissionRates", "input of the transmission rates"); 
    public Input<RealParameter> recoveryRatesInput = new Input<>("recoveryRates", "input of the recovery rates"); 
    public Input<RealParameter> recoveryRatesRatioInput = new Input<>("recoveryRatesRatio", "input of the recovery rates"); 
    public Input<RealParameter> waningImmunityRatesInput = new Input<>("waningImmunityRates", "input of the waning immunity rates"); 
    public Input<RealParameter> populationSizeInput = new Input<>("populationSize", "number of susceptibles at the beginning in each state");     
    public Input<RealParameter> initialIntroductionInput = new Input<>("initialIntroduction", "how many individuals existed at the root");
    public Input<RealParameter> samplingRatesInput = new Input<>("samplingRates", "input of the sampling rates");     
    public Input<RealParameter> migrationRatesInput = new Input<>("migrationRates", "input of the transmission rates");     
    public Input<BooleanParameter> MigrationIndicatorInput = new Input<>("migrationRatesIndicator", "indicates if a migration rate is zero or non zero");
    
    // specifies which states as SIR and some as constant
    public Input<BooleanParameter> stateIsConstantInput = new Input<>("stateIsConstant", "indicates if a state has constant rates");     

    // specify some other things   
    public Input<IntegerParameter> timeStepsInput = new Input<>("timeSteps", "number of time steps"); 
    public Input<RealParameter> initialNumberOfInfectedInput = new Input<>("initialNumberOfInfected", "initial number of infected for SIR simulation"); 
   
    
    private int states;
    public boolean hasChangingStates;
    
    private ArrayList<Double> time  = new ArrayList<>();
    private ArrayList<Double[]> birthRates  = new ArrayList<>();
    private ArrayList<Double[][]> migrationRates  = new ArrayList<>();
    private ArrayList<Double[]> infected  = new ArrayList<>();
    
    int[] stateMapping;
    int[] minStateMapping;
    
    
    private boolean areAllStatesConstant = false;
    private boolean hasSamplingRates = false;
    
	@Override
	public void initAndValidate() {	
		// check if all states are constant rate or if an SIR type model is to use
		if (constantCoalescentRatesInput.get() != null){
			areAllStatesConstant = true;
			states = constantCoalescentRatesInput.get().getDimension();
		}else{
			if (transmissionRatesInput.get().getDimension()==1)
				states = initialIntroductionInput.get().getDimension();
			else
				states = transmissionRatesInput.get().getDimension();
			
			if (samplingRatesInput.get() != null)
				hasSamplingRates = true;
		}
	}
	
	protected boolean update(double treeHeight, int[] stateMapping, int[] MinStateMapping){
		// get the state mapping from the parent class
		this.stateMapping = stateMapping;
		this.minStateMapping = MinStateMapping;
		
		// check if dynamics need to be calculated
		if(areAllStatesConstant){
			setConstantRates(treeHeight);
			return true;
		}else{
			return calculateSIRtrajectories(treeHeight);	
		}
	}
	
	/**
	 * Function builds the lists with the constant rates
	 * @param treeHeight
	 */
	private void setConstantRates(double treeHeight) {
		// set back array lists
		time = new ArrayList<>();
		birthRates = new ArrayList<>();
		migrationRates = new ArrayList<>();
		infected = new ArrayList<>();
		
		// build the time arrayList
		time.add(0.0);
		time.add(treeHeight/2);
		time.add(treeHeight);
		
		// get from the coalescent rates the birth rates 
		// (the divided by 2 is there to convert coalescent to birth rates)
		Double[] bRates = new Double[states];
		for (int i = 0; i < states; i++) 
			bRates[i] = constantCoalescentRatesInput.get().getArrayValue(i)/2;
		
		// set the number of infected to 1 in each state, so the birth rates translate
		// to the pairwise coalescent rates using a factor of 2
		Double[] inf = new Double[states];
		for (int i = 0; i < states; i++) 
			inf[i] = 1.0;
		
		// array for the migration rates
		Double[][] mRates = new Double[states][states];
		
		
		// Set the migration rates in the case of one global rate
		if (constantMigrationRatesInput.get().getDimension() == 1){
			for (int i = 0; i < states; i++){
				for (int j = 0; j < states; j++){
					if (i == j){
						mRates[i][j] = Double.POSITIVE_INFINITY;
					}else{
						mRates[i][j] = constantMigrationRatesInput.get().getValue();
					}
				}
			}		
		// Set the migration rates in the case of symmetric migration between states
		}else if (constantMigrationRatesInput.get().getDimension() == states*(states-1)/2){
			int c = 0; // counter for the migration rates
			for (int i = 0; i < states; i++){
				for (int j = 0; j < states; j++){
					if (i == j){
						mRates[i][j] = Double.POSITIVE_INFINITY;
					}else{
						if (constantMigrationIndicatorInput.get() != null){
							if (constantMigrationIndicatorInput.get().getArrayValue(c) == 0){
								mRates[i][j] = 0.0;
								c++;
							}else{
								mRates[i][j] = constantMigrationRatesInput.get().getArrayValue(c);
								c++;
							}						
						}else{
							mRates[i][j] = constantMigrationRatesInput.get().getValue(c);
							c++;
						}
					}
				}
			}	
		// Set the migration rates in the case of asymmetric migration
		}else if (constantMigrationRatesInput.get().getDimension() == states*(states-1)){
			int c = 0; // counter for the migration rates
			for (int i = 0; i < states; i++){
				for (int j = 0; j < states; j++){
					if (i == j){
						mRates[i][j] = Double.POSITIVE_INFINITY;
					}else{
						if (constantMigrationIndicatorInput.get() != null){
							if (constantMigrationIndicatorInput.get().getArrayValue(c) == 0){
								mRates[i][j] = 0.0;
								c++;
							}else{
								mRates[i][j] = constantMigrationRatesInput.get().getArrayValue(c);
								c++;
							}						
						}else{
							mRates[i][j] = constantMigrationRatesInput.get().getValue(c);
							c++;
						}
					}
				}
			}	
		}else {
			System.err.println("The dimension of the migration rates isn't chosen correctly choose either:\n"+ 
					"-dimension 1 for one global migration rate\n" + 
					"-dimension states*(state-1)/2 for symmetric migration\n" +
					"-dimension states*(state-1) for asymmetric migration\n");
		}	
		
		// build the birth rates, migration rates and the number of infected lists
		birthRates.add(bRates);birthRates.add(bRates);
		migrationRates.add(mRates);migrationRates.add(mRates);
		infected.add(inf);infected.add(inf);	
	}
	
	double[][] stateTime;
	double[][] stateBirth;
	double[][] stateInfected;

	private boolean calculateSIRtrajectories(double treeHeight) {
		// initialize the time arrays for each state, i.e.
		stateTime = new double[states][timeStepsInput.get().getValue()*15 + 1];		
		// initialize the birth rates in each state, i.e.
		stateBirth = new double[states][timeStepsInput.get().getValue()*15];
		// initialize the number of infected in each state, i.e.
		stateInfected = new double[states][timeStepsInput.get().getValue()*15];
		
		
		
		boolean isValid = false;
		for (int i = 0; i <  initialIntroductionInput.get().getDimension(); i++){
			if ( initialIntroductionInput.get().getArrayValue(i)>=treeHeight)
				isValid = true;
			if ( initialIntroductionInput.get().getArrayValue(i)>=treeHeight*2)
				isValid = false;
		}
		
		if (!isValid)
			return false;
		
		// initialize all rates
		double[][] mig_rates = getMigRates();
		double[] trans_rates = getTransRates();
		double[] rec_rates = getRecRates(trans_rates);
		double[] wan_rates = getWanRates();
		
        FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(treeHeight/(timeStepsInput.get().getValue()*5));
		// integrate until half the interval
		
		
		// set the initial number of infected and the initial time
		for (int i = 0; i < states; i++){
			FirstOrderDifferentialEquations ode = new SIR(trans_rates[i], rec_rates[i], wan_rates[i]);
			// initialize the things
			stateTime[i][0] = initialIntroductionInput.get().getArrayValue(stateMapping[i]);
			
			double[] p = new double[3];	
			double inititalInfected;
			
			if (initialNumberOfInfectedInput.get() != null){
				p[0] = populationSizeInput.get().getArrayValue(stateMapping[i]) - initialNumberOfInfectedInput.get().getValue();
				p[1] = initialNumberOfInfectedInput.get().getValue();
			}else{
				p[0] = populationSizeInput.get().getArrayValue(stateMapping[i]) * 0.9999;
				p[1] = populationSizeInput.get().getArrayValue(stateMapping[i]) * 0.0001;
			}			
			p[2] = 0.0;
			inititalInfected = p[1];
			double integrationStepSize = treeHeight/(timeStepsInput.get().getValue()*5);
			
			for (int j = 0; j < stateInfected[i].length; j++){
				if (p[1] >= inititalInfected){
					double[] p_for_ode = new double[3];	
					integrator.integrate(ode, 0, p, integrationStepSize, p_for_ode);						
					p = p_for_ode;				
					stateInfected[i][j] = p[1];
					stateBirth[i][j] = p[0]*p[1]/populationSizeInput.get().getArrayValue(stateMapping[i]) * trans_rates[i];
					stateTime[i][j+1] = stateTime[i][j] - integrationStepSize;
				}else{
					stateInfected[i][j] = 0.0;
					stateBirth[i][j] = 0.0;
					stateTime[i][j+1] = stateTime[i][j] - integrationStepSize;					
				}
			}			
		}
		
		
		// build the time list
		time = new ArrayList<>();
		double currTime = 0.0;
		for (int i = 0; i <= timeStepsInput.get().getValue()+1; i++){
			time.add(currTime);
			currTime += treeHeight/timeStepsInput.get().getValue();
		}
		birthRates = new ArrayList<>();
		infected = new ArrayList<>();

		// interpolate the above lists
		interpolate(stateInfected, stateBirth, stateTime);
		
		// reverse the ArraLists to go backward in time;
		reverseLists();
		buildMigrationRateList(mig_rates);
		updateConstantStates(mig_rates);
		return true;
	}


	
	private void interpolate(double[][] stateInfected, double[][] stateBirth, double[][] stateTime) {
		int[] start1 = new int[states];
		int[] start2 = new int[states];
		
		// very cheap way to interpolate
		for (int i = 0; i < states; i++){
			int c1 = time.size()-1;
			int c2 = 0;		
			
			// find the "initial" point
			if (stateTime[i][c2] > time.get(c1)){ // the state is "active" before the tree starts
				// the "target" time is the center of the first interval
				double targetTime = (time.get(c1-1) + time.get(c1))/2;
				double currentDifference = Math.abs(targetTime - stateTime[i][c2]);
				c2++;
				double nextDifference = Math.abs(targetTime - stateTime[i][c2]);
				
				while (nextDifference < currentDifference){
					currentDifference  = nextDifference;
					c2++;
					nextDifference = Math.abs(targetTime - stateTime[i][c2]);					
				}	
				c2--;	
				c1--;
				
			}else{// the state is inactive at the beginning and becomes active only later
				// find the first interval where the state is active
				while (time.get(c1) > stateTime[i][c2]) 
					c1--;
				
				double targetTime = (time.get(c1) + time.get(c1+1))/2;
				
				// find c2 such that the difference between stateTime[i][c2] and targetTime is minimal
				double currentDifference = Math.abs(targetTime - stateTime[i][c2]);
				c2++;
				double nextDifference = Math.abs(targetTime - stateTime[i][c2]);
				
				while (nextDifference < currentDifference){
					currentDifference  = nextDifference;
					c2++;
					nextDifference = Math.abs(targetTime - stateTime[i][c2]);					
				}	
				
			}
			start1[i] = c1;
			start2[i] = c2;
		}
	
		
		// build the lists of birth rates and the number of infected
		for (int i = time.size()-2; i >= 0 ; i--){
			Double[] bRates = new Double[states];
			Double[] nrInf = new Double[states];
			for (int j = 0; j < states; j++){
				if(start1[j] >= i){
					bRates[j] = stateBirth[j][start2[j]];
					nrInf[j] = stateInfected[j][start2[j]];
					start2[j] += 5;
				}else{
					bRates[j] = 0.0;
					nrInf[j] = 0.0;
				}
				
			}
//			System.out.println("time = " + time.get(i) + " number of infected = " +Arrays.toString(nrInf));
			birthRates.add(bRates);
			infected.add(nrInf);
		}
	}

	/*
	 * Update all the rates given the Input values and an SIR type model
	 */
	private double[][] getMigRates(){
		double[][] mig_rates = new double[states][states];
		// fill the rates arrays: migration rates
		int c = 0;
		if (migrationRatesInput.get().getDimension()==1){
			for (int i = 0; i < states; i++){
				for (int j = 0; j < states; j++){
					if(i!=j){
						mig_rates[i][j] = migrationRatesInput.get().getValue();
					}else{
						mig_rates[i][j] = 0.0;
					}
				}
			}
		}else{		
			for (int i = 0; i < states; i++){
				for (int j = 0; j < states; j++){
					if(i!=j){
						if (MigrationIndicatorInput.get()!=null){
							if (MigrationIndicatorInput.get().getArrayValue(c)==1.0)
								mig_rates[i][j] = migrationRatesInput.get().getArrayValue(c);
							else
								mig_rates[i][j] = 0.0;
						}else{
							mig_rates[i][j] = migrationRatesInput.get().getArrayValue(c);
						}
						c++;
					}else{
						mig_rates[i][j] = 0.0;
					}
				}
			}
		}
		return mig_rates;
	}
	
	private double[] getTransRates(){
		double[] trans_rates = new double[states];
		for (int i = 0; i < states; i++){
			if(transmissionRatesInput.get().getDimension()==1){
				trans_rates[i] = transmissionRatesInput.get().getValue();
			}else{
				trans_rates[minStateMapping[i]] = transmissionRatesInput.get().getArrayValue(stateMapping[i]);
			}
		}
		return trans_rates;
	}
	
	private double[] getRecRates(double[] trans_rates){
		double[] rec_rates = new double[states];
		if (recoveryRatesInput.get()!=null)
			if(recoveryRatesInput.get().getDimension()==1)
				for (int i = 0; i < states; i++)
					rec_rates[i] = recoveryRatesInput.get().getValue();
			else
				for (int i = 0; i < states; i++)
					rec_rates[minStateMapping[i]] = recoveryRatesInput.get().getArrayValue(stateMapping[i]);
		else
			for (int i = 0; i < states; i++)
				rec_rates[i] = 0.0;
		
		if (recoveryRatesRatioInput.get()!=null){
			for (int i = 0; i < states; i++)
				if (recoveryRatesRatioInput.get().getDimension()==1)					
					rec_rates[i] = recoveryRatesRatioInput.get().getValue()*trans_rates[i];
				else
					rec_rates[minStateMapping[i]] = recoveryRatesRatioInput.get().getArrayValue(i)*trans_rates[stateMapping[i]];
		}

		return rec_rates;
	}
	
	private double[] getWanRates(){
		double[] wan_rates = new double[states];
		if (waningImmunityRatesInput.get()!=null)
			for (int i = 0; i < states; i++)
				wan_rates[minStateMapping[i]] = waningImmunityRatesInput.get().getArrayValue(stateMapping[i]);
		else
			for (int i = 0; i < states; i++)
				wan_rates[i] = 0.0;

		return wan_rates;
	}	
	
	private void updateConstantStates(double[][] mig_rates){
		if(stateIsConstantInput.get() == null)
			return;
					
		// check if there are states with constant rates
		boolean hasConstState = false;
		for (int i = 0; i < states; i++){
			// Boolean parameter returns double
			if (stateIsConstantInput.get().getArrayValue(i) == 1)
				hasConstState = true;
		}
		
		// if there is nothing to do return
		if (!hasConstState)
			return;
		
		// if there are constant states, get the mean number of infected for that state
		double[] nrInfected = new double[states];
		for (int i = 0; i < states; i++){
			for (int j = 0; j < infected.size(); j++)
				nrInfected[i] += infected.get(j)[i];
			// take the mean number of infected
			nrInfected[i] /= infected.size();
		}
		
		// newly set the birth rates and the number of infected for the constant states
		for (int i = 0; i < states; i++){
			// Boolean parameter returns double
			if (stateIsConstantInput.get().getArrayValue(i) == 1){
				for (int j = 0; j < infected.size(); j++){
					infected.get(j)[i] = nrInfected[i];
					birthRates.get(j)[i] = transmissionRatesInput.get().getArrayValue(i);
				}				
			}
		}
		
		// update the migration rates (change in the future as it seems the wrong place to do that)
		if (states > 1){
			for (int i = 0; i < states; i++){
				for (int j = 0; j < states; j++){
					for (int k = 0; k < migrationRates.size(); k++){
						if (i!=j){
							if (infected.get(k)[j]>0){		
								migrationRates.get(k)[j][i] = infected.get(k)[i] * mig_rates[i][j] / infected.get(k)[j];
							}else{
								migrationRates.get(k)[j][i] = Double.POSITIVE_INFINITY;
							}
						}else{
							migrationRates.get(k)[i][j] = Double.POSITIVE_INFINITY;
						}
					}
				}
			}
		}
	}

	/**
	 * Function reverses all the ArraysList such that they are backwards in time for the coalescent
	 */
	private void reverseLists(){
//		ArrayList<Double> newTime = new ArrayList<>();
		ArrayList<Double[]> newBirthRates = new ArrayList<>();
		ArrayList<Double[]> newNumberOfInfected = new ArrayList<>();
		
		// get the list in reverse
//		for (int i = time.size()-1; i>=0; i--){
//			newTime.add(time.get(i));
//		}
		for (int i = birthRates.size()-1; i>=0; i--){
			newBirthRates.add(birthRates.get(i));
			newNumberOfInfected.add(infected.get(i));
		}

		// set new to old
//		time = new ArrayList<>(newTime);
		birthRates = new ArrayList<>(newBirthRates);
		infected = new ArrayList<>(newNumberOfInfected);
	}
	
	private void buildMigrationRateList(double[][] mig_rates){
		migrationRates = new ArrayList<>();
		for (int i = 0; i < birthRates.size(); i++){
			Double[][] newMigRates = new Double[states][states];
			
			if (states>1){
				for (int j = 0; j < states; j++){
					for (int k = 0; k < states; k++){
						if (j == k){
							newMigRates[j][k] = Double.POSITIVE_INFINITY;						
						}else{
							if (infected.get(i)[j] > 0)
								newMigRates[j][k] = mig_rates[j][k] * infected.get(i)[k]/infected.get(i)[j];
							else
								newMigRates[j][k] = Double.POSITIVE_INFINITY;
						}
					}
				}
			}else{
				newMigRates = new Double[1][1];
				newMigRates[0][0] = Double.POSITIVE_INFINITY;
			}
			

			
			migrationRates.add(newMigRates);
//			System.out.println();
//			System.out.print("[");
//			for (int a = 0; a < states; a++){
//				for (int b = 0; b < states-1; b++){
//					System.out.print(migrationRates.get(0)[a][b] + ",");			
//				}
//				System.out.print(migrationRates.get(0)[a][states-1] + ";");			
//			}
//			System.out.print("]\n");

		}
	}

	/**
	 * Function that return 2 if all states are constant and the input number of steps elsewise 
	 * @return
	 */
	public int getNrSteps(){
		if (areAllStatesConstant)
			return 2;
		else
			return timeStepsInput.get().getValue();
	}
	
	
	public ArrayList<Double> getTime(){
		return time;		
	}
	
	public ArrayList<Double[]> getBirthRates(){
		return birthRates;		
	}
	
	public ArrayList<Double[]> getNumberOfInfectedRates(){
		return infected;		
	}
	
	public ArrayList<Double[][]> getMigrationRates(){
		return migrationRates;		
	}
			
	protected boolean hasSamplingRates(){
		return hasSamplingRates;
	}
	
	public int getStates(){
		return states;
	}
	
	@Override
	public void init(PrintStream out) {
		for (int i = 0 ; i < 4; i++){
			for (int j = 0; j < birthRates.size(); j++)
				out.print("Ne" + i + "." + j + "\t");
		}
//		for (int i = 0 ; i < 6; i++){
//			for (int j = 0; j < birthRates.size(); j++)
//				out.print("Ne" + i + "." + j + "\t");
//		}



	}

	@Override
	public void log(int sample, PrintStream out) {
		for (int i = 0 ; i < states; i++){
			for (int j = 0; j < birthRates.size(); j++)
				if ((infected.get(j)[i]*infected.get(j)[i])/birthRates.get(j)[i] >= 0.0)
					out.print((infected.get(j)[i]*infected.get(j)[i])/birthRates.get(j)[i] + "\t");
				else
					out.print("0.0\t");
		}
		for (int i = 0 ; i < states; i++){
			for (int j = 0; j < birthRates.size(); j++)
				if ((infected.get(j)[i]*infected.get(j)[i]) >= 0.0)
					out.print((infected.get(j)[i]*infected.get(j)[i])/birthRates.get(j)[i] + "\t");
				else
					out.print("0.0\t");
		}

//		for (int i = 1 ; i < states; i++){
//			for (int j = 0; j < birthRates.size(); j++)
////				if ((infected.get(j)[i]*infected.get(j)[i])/birthRates.get(j)[i] >= 0.0)
//					out.print(migrationRates.get(j)[1][0] + "\t");
////				else
////					out.print("5.0\t");
//		}
//		for (int i = 1 ; i < states; i++){
//			for (int j = 0; j < birthRates.size(); j++)
////				if ((infected.get(j)[i]*infected.get(j)[i])/birthRates.get(j)[i] >= 0.0)
//					out.print(migrationRates.get(j)[0][1] + "\t");
////				else
////					out.print("5.0\t");
//		}
//
//
//		for (int i = 2 ; i < 3; i++){
//			for (int j = 0; j < birthRates.size(); j++)
////				if ((infected.get(j)[1]*infected.get(j)[1])/birthRates.get(j)[1] >= 0.0)
//					out.print(infected.get(j)[0]/infected.get(j)[1] + "\t");
////				else
////					out.print("0.0\t");
//		}
//		
//		for (int i = 1 ; i < states; i++){
//			for (int j = 0; j < birthRates.size(); j++)
//				if (infected.get(j)[i] == 0.0)
//					out.print("0.0\t");
//				else
//					out.print("1.0\t");
//		}
//		
//			for (int j = 0; j < birthRates.size(); j++)
//				out.print(infected.get(j)[0] + "\t");
//			for (int j = 0; j < birthRates.size(); j++)
//				out.print(infected.get(j)[1] + "\t");
		



	}

	@Override
	public void close(PrintStream out) {
		out.close();
	}
	
	@Override
    protected boolean requiresRecalculation() {
    	return true;
    }   

}
