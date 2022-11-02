package mascot.mapped;



import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.IntervalList;
import beast.base.evolution.tree.IntervalType;

/**
 * Extracts the intervals from a beast.tree. *
 * @author Nicola Felix Müller
 */
@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class MappedTreeIntervals extends CalculationNode {
    public Input<Tree> treeInput = new Input<Tree>("tree", "tree for which to calculate the intervals", Validate.REQUIRED);
    public Input<Double> maxIntervalLengthInput = new Input<Double>("maxIntervalLength",
    		"the maximal length of an interval in absolute time before its split into different storing points");
    public Input<Double> maxRelativeIntervalLengthInput = new Input<Double>("maxRelativeIntervalLength",
    		"the maximal length of an interval relative to the tree height before its split into different storing points", Validate.XOR, maxIntervalLengthInput);
   
    
    // keeps track of the different types of events needed for state clocks to work
    public enum Events{Coalescent, Sampling, Storing};
    
    private Double maxIntervalLength;
    
    private Tree tree;
    /**
     * The widths of the intervals.
     */
    public double[] intervals;
   
    /**
     * The number of uncoalesced lineages within a particular interval.
     */
    protected Type[] types;
    
    protected boolean intervalsKnown = false;

    /**
     * =======================================================================
     * =======================================================================
     * =======================================================================
     * =======================================================================
     */
    
    public MappedTreeIntervals() {
    }

    public MappedTreeIntervals(Tree tree) {
    	calculateIntervals();
    }

    @Override
    public void initAndValidate() {
    	calculateIntervals();
    }

    /**
     * CalculationNode methods *
     */
    @Override
    protected boolean requiresRecalculation() {
        // we only get here if the tree is dirty, which is a StateNode
        // since the StateNode can only become dirty through an operation,
        // we need to recalculate tree intervals
        intervalsKnown = false;
        return true;
    }
    
    public int getSampleCount() {
        // Assumes a binary tree!
        return tree.getInternalNodeCount();
    }
    

    /**
     * get number of intervals
     */
    public int getIntervalCount() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals.length;
    }

    /**
     * Gets an interval.
     */
    public double getInterval(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < 0 || i >= intervals.length){
        	throw new IllegalArgumentException();
        }
        return intervals[i];
    }
    
    /**
     * Defensive implementation creates copy
     *
     * @return
     */
    public double[] getIntervals(double[] inters) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (inters == null) inters = new double[intervals.length];
        System.arraycopy(intervals, 0, inters, 0, intervals.length);
        return inters;
    }
    
   
    
    public int getLineagesAdded(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervals.length) throw new IllegalArgumentException();
        
        
        return types[i].getNodeNr();
    }
    
    public int getLineagesRemoved(int index, int index2) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (index >= intervals.length) throw new IllegalArgumentException();
        
        if (index2==0)
        	return tree.getNode(types[index].getNodeNr()).getLeft().getNr();
        else
        	return tree.getNode(types[index].getNodeNr()).getRight().getNr();
        
    }

    /**
     * Returns the type of interval observed.
     */
    public Events getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervals.length) throw new IllegalArgumentException();
        return types[i].getEvent();
    }


    /**
     * Recalculates all the intervals for the given beast.tree.
     */
    @SuppressWarnings("unchecked")
	public void calculateIntervals() {
		tree = treeInput.get();
		
		// calculate the maximimal length of an interval before it has to be split apart
		if (maxRelativeIntervalLengthInput.get()!=null)
			maxIntervalLength = maxRelativeIntervalLengthInput.get()*tree.getRoot().getHeight();
		else
			maxIntervalLength = maxIntervalLengthInput.get();
		
		List<Double> eventTimes = new ArrayList<>();
		List<Type> eventTypes = new ArrayList<>();
		List<Node> nodes = new ArrayList<>();

        collectTimes(tree, eventTimes, eventTypes, maxIntervalLength, nodes);

        // calculate the sorted tree intervals
        intervals = new double[eventTimes.size()];
        types = new Type[eventTimes.size()];
        
        //get the sorted event times
        List<Double> sortedTimes = new ArrayList<>(eventTimes);
        Collections.sort(sortedTimes);
        
        double currentTime = 0.0;

        for (int i = 0; i < sortedTimes.size(); i++){
        	// check if there are more than one coalescent events with the same event time, if yes
        	// ensure that how the nodes are connected is actually reflected by the
        	int index = eventTimes.indexOf(sortedTimes.get(i));
        	int lastIndex = eventTimes.lastIndexOf(sortedTimes.get(i));
        	
        	if (index!=lastIndex) {
        		int j=i;
        		double currentEventTime = sortedTimes.get(i);
        		List<Type> coalTypes = new ArrayList<>();
        		List<Node> coalNodes = new ArrayList<>();

        		while (index!=-1) {
        			if (eventTypes.get(index).getEvent()==Events.Coalescent) {
        				coalTypes.add(eventTypes.get(index));
        				coalNodes.add(nodes.get(index));
        	        	eventTimes.set(index, Double.NaN);
        			}else {
        	        	types[j] = eventTypes.get(index);
        	        	intervals[j] = currentEventTime - currentTime;
        	        	currentTime = currentEventTime;
        	        	eventTimes.set(index, Double.NaN);
        	        	j++;
        			}
        			index = eventTimes.indexOf(currentEventTime);
        		}
        		
        		// if there is only one coalescent event
        		if (coalTypes.size()==1) {
    	        	types[j] = coalTypes.get(0);
    	        	intervals[j] = currentEventTime - currentTime;
    	        	j++;
        		} else if (coalTypes.size()==2) { // if there are two coalescent events at the same time
        			if (coalNodes.get(0).getParent()==coalNodes.get(1)){ 
	    	        	types[j] = coalTypes.get(0);
	    	        	intervals[j] = currentEventTime - currentTime;
	    	        	j++;
	    	        	types[j] = coalTypes.get(1);
	    	        	intervals[j] = currentEventTime - currentTime;
	    	        	j++;
	    	        	currentTime = currentEventTime;
        			}
        			else {
	    	        	types[j] = coalTypes.get(1);
	    	        	intervals[j] = currentEventTime - currentTime;
	    	        	j++;
	    	        	types[j] = coalTypes.get(0);
	    	        	intervals[j] = currentEventTime - currentTime;
	    	        	j++;
	    	        	currentTime = currentEventTime;
        			}

        		}else {
	        		// for any other case, explicitly check the orderings
        			// TODO actually check the ordering
	        		for (int k = 0; k<coalTypes.size(); k++) {
	    	        	types[j] = coalTypes.get(k);
	    	        	intervals[j] = currentEventTime - currentTime;
	    	        	j++;
	    	        	currentTime = currentEventTime;
	        		}
        		}
        		i=j-1;
        		
        	}else {
	        	// get the corresponding indices
	        	types[i] = eventTypes.get(index);
	        	intervals[i] = sortedTimes.get(i) - currentTime;
	        	currentTime = sortedTimes.get(i);
	        	eventTimes.set(index, Double.NaN);
        	}
        }
        
        intervalsKnown = true;
    }


    /**
     * extract coalescent times and tip information into array times from beast.tree.
     *
     * @param tree        the beast.tree
     * @param times       the times of the nodes in the beast.tree
     * @param childCounts the number of children of each node
     */
    protected void collectTimes(Tree tree, List<Double> times, List<Type> types, double maxLength, List<Node> nodesList) {
        Node[] nodes = tree.getNodesAsArray();
                
        for (int i = 0; i < nodes.length; i++) {
            Node node = nodes[i];
            // check the length of the node.
            if (!node.isRoot()){
            	double length = node.getParent().getHeight() - node.getHeight();
            	if (length>maxLength){
            		// add the current node
            		times.add(node.getHeight());
            		
            		if (node.isLeaf()) {
            			types.add(new Type(node.getNr(), Events.Sampling));
                		nodesList.add(null);
            		}
            		else {
            			types.add(new Type(node.getNr(), Events.Coalescent));
                		nodesList.add(node);
            		}

            		// split the node into several intervals
            		int nrIntervals = (int) (length/maxLength) + 1;
            		for (int j = 1; j < nrIntervals-1; j++){
                		times.add(node.getHeight() + j* length/nrIntervals);
                		types.add(new Type(node.getNr(), Events.Storing));
                		nodesList.add(null);
            		}
            	}else{
                	times.add(node.getHeight());
                	nodesList.add(node);
            		if (node.isLeaf()) types.add(new Type(node.getNr(), Events.Sampling));
            		else types.add(new Type(node.getNr(), Events.Coalescent));
            	}
            }else{
            	times.add(node.getHeight());
        		types.add(new Type(node.getNr(), Events.Coalescent));
        		nodesList.add(node);
            }            
        }
    }
    
    private class Type {
    	Type(int nodeNr, Events event){
    		this.nodeNr = nodeNr;
    		this.event = event;    				
    	}
    	
    	int getNodeNr(){
    		return nodeNr;
    	}
    	
    	Events getEvent(){
    		return event;
    	}   	
    	
    	int nodeNr;
    	Events event;
    }
}