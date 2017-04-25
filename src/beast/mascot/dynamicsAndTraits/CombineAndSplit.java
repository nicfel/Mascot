/*
* File ScaleOperator.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
package beast.mascot.dynamicsAndTraits;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.parameter.RealParameterList;
//import beast.core.parameter.BooleanParameter;
//import beast.core.parameter.IntegerParameter;
//import beast.evolution.tree.Node;
//import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
 * 
 * changed by Nicola Felix Mueller from the scale operator
 *
 */
@Description("Adds a fixed value or multiplies by a fixed value. Is always accepted")
public class CombineAndSplit extends Operator {
	
    public Input<IntegerParameter> stateMappingInput = new Input<>(
    		"stateMapping",
    		"how states are mapped"); 
    
    public Input<IntegerParameter> nrStatesInput = new Input<>(
    		"nrStates",
    		"how states are mapped"); 
 
    // List of 1 and 2 dimensional parameters that eventually have to be set to be the same
    final public Input<List<RealParameter>> parameter1D =
            new Input<>("parameter1D",
                    "parameters that have to be set equal and have only one dimension",
                    new ArrayList<>());
    final public Input<List<RealParameter>> parameter2D =
            new Input<>("parameter2D",
                    "parameters that have to be set equal and have only one dimension",
                    new ArrayList<>());

    
    ArrayList<Integer> currentStates;
    ArrayList<Integer> numberOfStates;
    
    int[] stateMap;
    
    int[] numberOfRealParametersInput;
    
    @Override
    public void initAndValidate() {
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	// calculate which states are active and how many sub states
    	// each state has
    	getCurrentStates();
    	
    	// get the input states
    	IntegerParameter nrStatesParam = nrStatesInput.get();
    	
    	// split
		if (Randomizer.nextBoolean()){
    		if (currentStates.size() == stateMappingInput.get().getDimension()){
    			return Double.NEGATIVE_INFINITY;
    		}
    		split();
			nrStatesParam.setValue(nrStatesParam.getValue()+1);
    	}
    	else{
    		if (currentStates.size()==1){
    			return Double.NEGATIVE_INFINITY;
    		}
    		combine();
			nrStatesParam.setValue(nrStatesParam.getValue()-1);
    	}
		// update the State mapping parameter
    	IntegerParameter stateMappingParam = stateMappingInput.get();
    	for (int i = 0; i < stateMappingParam.getDimension(); i++)
    		stateMappingParam.setValue(i, stateMap[i]);
    	
    	return 0.0;
    	
    }
    
    private void getCurrentStates(){
    	// build the state map
    	stateMap = new int[stateMappingInput.get().getDimension()];
    	for (int i = 0; i < stateMap.length; i++)
    		stateMap[i] = (int) stateMappingInput.get().getArrayValue(i);
    	   	
    	currentStates = new ArrayList<>();
    	numberOfStates = new ArrayList<>();
    	for (int i = 0; i < stateMappingInput.get().getDimension(); i++)
    		currentStates.add((int) stateMappingInput.get().getArrayValue(i));

    	// make a copy of the current states
    	ArrayList<Integer> allStates = new ArrayList<>(currentStates);

    	// sort the current States
    	Collections.sort(currentStates);
    	Collections.sort(allStates);
    	
    	// remove double entries from current states
    	for (int i = currentStates.size()-1; i > 0; i--)
    		if (currentStates.get(i) == currentStates.get(i-1))
    			currentStates.remove(i);
    		
    	// get the number of entries of all current state
    	numberOfStates.add(1);
    	int counter = 0;
    	for (int i = 1; i < allStates.size(); i++){
    		if (allStates.get(i) == allStates.get(i-1)){
    			numberOfStates.set(counter, numberOfStates.get(counter) + 1);
    		}else{
    			counter++;
    	    	numberOfStates.add(1);
    		}
    	}
    }
    
    private void split(){
    	// define which state to split
        int splitState = Randomizer.nextInt(currentStates.size());
        while(numberOfStates.get(splitState) == 1)
        	splitState = Randomizer.nextInt(numberOfStates.size());
        
        // initialize the array that knows how to split the states
        boolean[] targetState = new boolean[numberOfStates.get(splitState)];
       
        // inefficient way to split in exactly two parts
        int hasNrSplits = 0;
        do {
            hasNrSplits = 0;
            for (int i = 0; i < targetState.length; i++){
            	if (Randomizer.nextBoolean()){
            		targetState[i] = true;
            		hasNrSplits++;
            	}else{
            		targetState[i] = false;         		
            	}
            }
        	
        }while((targetState.length-1) < hasNrSplits || hasNrSplits == 0);
      	splitState(splitState, targetState);
    }
    
    private void combine(){
    	// take two random states which are to combine
        int State1 = Randomizer.nextInt(currentStates.size());
        int State2 = Randomizer.nextInt(currentStates.size());
        while (State1==State2)
        	State2 = Randomizer.nextInt(currentStates.size());     
        
        combineStates(currentStates.get(State1), currentStates.get(State2));
        combineParameters(currentStates.get(State1), currentStates.get(State2));
    }

 	private void combineStates(int val1, int val2){
    	
		int idx1 = currentStates.indexOf(val1);
		int idx2 = currentStates.indexOf(val2);		
   	
    	// Set the value of state val2 to val1 or vice versa (lower first)
    	if (idx1 < idx2){
    		for (int i = 0; i < stateMap.length; i++)
    			if (stateMap[i] == val2)
    				stateMap[i] = val1;

    	}else if(idx2 < idx1){
    		for (int i = 0; i < stateMap.length; i++)
    			if (stateMap[i] == val1)
    				stateMap[i] = val2;    	    		
   		
    	}else{
    		System.out.println("value 1 and value 2 are the same, can't combine the same state\n"+
    							"exit now:...");
    		System.exit(0);
    	}
  }
    
 	private void combineParameters(int val1, int val2){
 		if (parameter1D.get() != null){
 			int[] minStateMap = getMinimalStateMapping();
 			for (int i = 0; i < parameter1D.get().size(); i++){
 				RealParameter param = parameter1D.get().get(i);
 				double value1 = param.getArrayValue(val1);
 				double value2 = param.getArrayValue(val2);
 				
 				// take the mean of these values
 				double newValue = (value1+value2)/2;
 				
 				// set the new parameter values
 				for (int j = 0; j < minStateMap.length; j++){
 					if (val1 == stateMap[j])
 						param.setValue(j, newValue);
 					if (val2 == stateMap[j])
 						param.setValue(j, newValue);					
 				}
 			}
 		}
 	}
   
 	
 	public int[] getMinimalStateMapping(){
    	// make an ArrayList copy of the trait mapping
    	ArrayList<Integer> traitMapCopy = new ArrayList<Integer>();
    	for (int i = 0; i < stateMap.length; i++)
    		traitMapCopy.add(stateMap[i]);
    	

    	// sort the current states
    	Collections.sort(traitMapCopy);
    	
    	int[] minStateMapping = new int[stateMap.length];
    	for (int i = traitMapCopy.size()-2; i >= 0; i--)
    		if (traitMapCopy.get(i) == traitMapCopy.get(i+1))
    			traitMapCopy.remove(i);
    	
    	for (int i = 0; i < stateMap.length; i++)
    		for (int j = 0; j < traitMapCopy.size(); j++)
    			if (traitMapCopy.get(j) == stateMap[i])
    				minStateMapping[i] = j;    		
    	
    	
    	return minStateMapping;
    }

    public void splitState(int val, boolean[] targetState){
		if (numberOfStates.get(val) != targetState.length){
    		System.out.println("there has been an issue in the number of sub states and the number of target states"+
    							" when trying to split a state\nexit now:...");
    		System.exit(0);			
		}
		
		// Get all the subState of val
		int[] hasSubStates = new int[targetState.length];
		int counter = 0;
		for (int i = 0 ; i < stateMap.length; i++){
			if(stateMap[i] == currentStates.get(val)){
				hasSubStates[counter] = i;
				counter++;
			}
		}
		
		if (counter==0){
			System.out.println("Error in the sub states assignment");
			System.out.println("exit now...");
			System.exit(0);
		}
		
		// find the lowest integer for both splitting groups
		int min1 = Integer.MAX_VALUE, min2 = Integer.MAX_VALUE;		
    	int minusStates = 0;
		for (int i = 0; i < targetState.length; i++){
			if(targetState[i]){
				if (hasSubStates[i] < min1)
					min1 = hasSubStates[i];
			}else{
				minusStates++;
				if (hasSubStates[i] < min2)
					min2 = hasSubStates[i];
			}
		}
		
		if (min1 == Integer.MAX_VALUE || min2 == Integer.MAX_VALUE){
			System.out.println("cant see a split, likely all states are proposed to split\nexit now...");
			System.out.println(Arrays.toString(stateMap));
			System.out.println(currentStates);
			System.out.println(Arrays.toString(hasSubStates));
			System.out.println(val);
			System.exit(0);			
		}
		
		// update the trait mapping
		for (int i = 0; i < targetState.length; i++){
			if (targetState[i]){
				stateMap[hasSubStates[i]] = min1;
			}else{
				stateMap[hasSubStates[i]] = min2;
			}
		}
   }

} // class ScaleOperator
