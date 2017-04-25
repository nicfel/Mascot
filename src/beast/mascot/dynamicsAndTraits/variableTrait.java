package beast.mascot.dynamicsAndTraits;

import java.io.PrintStream;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.w3c.dom.Node;

import beast.core.BEASTObject;
import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.Loggable;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;

/**
 * 
 * @author Nicola Felix Mueller
 *
 */

@Description("A trait class where individual traits can be combined or splitted")
public class variableTrait extends CalculationNode {
	
    public Input<TraitSet> typeTraitInput = new Input<>(
            "typeTrait", "Type trait set.");     
    
    public Input<IntegerParameter> stateMappingInput = new Input<>(
    		"stateMapping",
    		"how states are mapped"); 
    
    public Input<IntegerParameter> nrStatesInput = new Input<>(
    		"nrStates",
    		"how states are mapped"); 
   
    double[] values;
    double minValue;
    int maxValue;
            
    // map several traits to one etc.
    int[] stateMapping;
    int[] minStateMapping;

    /**
     * Whether or not values are ALL numeric.
     */
    boolean numeric = true;
    
    @Override
    public void initAndValidate() {
    	// set the maximum value
    	maxValue = stateMappingInput.get().getDimension();
    	
    	// check validity of trait mapping
    	for (int i = 0; i < stateMappingInput.get().getDimension(); i++){
    		if (stateMappingInput.get().getArrayValue(i) != i){
    			System.err.println("initial trait map does not follow 0,...,n-1, please change");
    			System.exit(1);
    		}
    	}	
        
        // initialize the state mapping
        stateMapping = new int[maxValue];
        for (int i = 0; i < maxValue; i++){
        	stateMapping[i] = (int) stateMappingInput.get().getArrayValue(i);
        }
        
       
    } // initAndValidate
    
    public void update(){
    	stateMapping = new int[stateMapping.length];
    	for (int i = 0; i < stateMappingInput.get().getDimension(); i++)
    		stateMapping[i] = (int) stateMappingInput.get().getArrayValue(i);
    	
    	minStateMapping = getMinimalStateMapping();   	
    }
    
    public int[] getMinimalStateMapping(){
    	// make an ArrayList copy of the trait mapping
    	ArrayList<Integer> traitMapCopy = new ArrayList<Integer>();
    	for (int i = 0; i < stateMapping.length; i++)
    		traitMapCopy.add(stateMapping[i]);
    	


    	// sort the current states
    	Collections.sort(traitMapCopy);
    	
  	
    	    	
    	minStateMapping = new int[stateMapping.length];
    	for (int i = traitMapCopy.size()-2; i >= 0; i--)
    		if (traitMapCopy.get(i) == traitMapCopy.get(i+1))
    			traitMapCopy.remove(i);
    	
    	for (int i = 0; i < stateMapping.length; i++)
    		for (int j = 0; j < traitMapCopy.size(); j++)
    			if (traitMapCopy.get(j) == stateMapping[i])
    				minStateMapping[i] = j;    		
    	
    	return minStateMapping;
    }

    public int[] getMinimalStateMapping(boolean update){
    	if(update)
    		update();
    	// make an ArrayList copy of the trait mapping
    	ArrayList<Integer> traitMapCopy = new ArrayList<Integer>();
    	for (int i = 0; i < stateMapping.length; i++)
    		traitMapCopy.add(stateMapping[i]);
    	


    	// sort the current states
    	Collections.sort(traitMapCopy);
    	
  	
    	    	
    	minStateMapping = new int[stateMapping.length];
    	for (int i = traitMapCopy.size()-2; i >= 0; i--)
    		if (traitMapCopy.get(i) == traitMapCopy.get(i+1))
    			traitMapCopy.remove(i);
    	
    	for (int i = 0; i < stateMapping.length; i++)
    		for (int j = 0; j < traitMapCopy.size(); j++)
    			if (traitMapCopy.get(j) == stateMapping[i])
    				minStateMapping[i] = j;    		
    	
    	return minStateMapping;
    }

    
    public int getSamplingState(String taxonName){
		int sampleState = (int) typeTraitInput.get().getValue(taxonName);
		return minStateMapping[sampleState];
    }
    
	@Override
	protected void store() {
	}

	@Override
	public void restore() {
	}

	public int[] getStateMapping() {
		return stateMapping;
	}

	public int getNrStates() {
		return nrStatesInput.get().getValue();
	}

    
} // class TraitSet
