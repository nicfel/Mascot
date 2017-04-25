package beast.mascot.dynamicsAndTraits;

import java.util.ArrayList;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;

/**
 * @author Nicola Felix Mueller 
 */
@Description("class that connects the rates and traits with the actual Likelihood calculations")
public class RatesAndTraits extends CalculationNode {	
    public Input<variableTrait> variableTraitInput = new Input<>("variableTrait", "traits that can change");
    public Input<Rates> ratesInput = new Input<>("rates", "traits that can change");
   
    public boolean hasChangingStates;
    
    private int[] stateMapping;
    private int[] minStateMapping;
    
    private boolean hasVariableTraits = false;    
    
	@Override
	public void initAndValidate() {
		// check if there are variable traits given as input
		if(variableTraitInput.get() != null) hasVariableTraits=true;
		
		// Initialize the state maps
		stateMapping = new int[ratesInput.get().getStates()];
		minStateMapping = new int[ratesInput.get().getStates()];
		for (int i = 0; i < ratesInput.get().getStates(); i++){
			stateMapping[i] = i;
			minStateMapping[i] = i;
		}

	}
	
	public boolean update(double treeHeight) {	
		// If variable traits are given as Input, use them
		if (hasVariableTraits){
			variableTraitInput.get().update();
			stateMapping = variableTraitInput.get().getStateMapping();
			minStateMapping = variableTraitInput.get().getMinimalStateMapping();
		}
		return ratesInput.get().update(treeHeight, stateMapping, minStateMapping);
	}

	public int getNrSteps(){
		return ratesInput.get().getNrSteps();
	}
	
	public ArrayList<Double> getTime(){
		return ratesInput.get().getTime();		
	}
	
	public ArrayList<Double[]> getBirthRates(){
		return ratesInput.get().getBirthRates();		
	}
	
	public ArrayList<Double[]> getNumberOfInfectedRates(){
		return ratesInput.get().getNumberOfInfectedRates();		
	}
	
	public ArrayList<Double[][]> getMigrationRates(){
		return ratesInput.get().getMigrationRates();		
	}		
		
	public int getSampleState(String taxonName){
		return variableTraitInput.get().getSamplingState(taxonName);	
	}

	public boolean hasSamplingRates(){
		return ratesInput.get().hasSamplingRates();
	}
	
	public double getSamplingRate(int state){
		return ratesInput.get().samplingRatesInput.get().getArrayValue(state);
	}
	
	public int getStates(){
		return ratesInput.get().getStates();
	}
		
	@Override
    protected boolean requiresRecalculation() {
    	return true;
    }   

	public boolean hasVariableTraits(){
		return hasVariableTraits;
	}

	
}
