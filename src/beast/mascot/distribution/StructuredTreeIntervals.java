package beast.mascot.distribution;



import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.IntervalList;
import beast.evolution.tree.coalescent.IntervalType;
import beast.util.HeapSort;


/*
 * TreeIntervals.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
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

/**
 * Extracts the intervals from a beast.tree.
 * DAR: modified this class to work with structured coalecent models
 * changed by Nicola Felix Mueller to make it compatible with the StructCoalDensity
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeIntervals.java,v 1.9 2005/05/24 20:25:56 rambaut Exp $
 */
@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class StructuredTreeIntervals extends CalculationNode implements IntervalList {
    public Input<Tree> treeInput = new Input<Tree>("tree", "tree for which to calculate the intervals", Validate.REQUIRED);

    private Tree tree;
    /**
     * The widths of the intervals.
     */
    protected double[] intervals;
    protected double[] storedIntervals;
    
   
    /**
     * If interval was changed 
     */
    protected boolean[] intervalIsDirty;

    /**
     * The number of uncoalesced lineages within a particular interval.
     */
    protected int[] lineageCounts;
    protected int[] storedLineageCounts;

    /**
     * The lineages in each interval (stored by node ref).
     */
    protected List<Node>[] lineagesAdded;
    protected List<Node>[] lineagesRemoved;
    
    // Added these so can restore when needed for structured models
    protected List<Node>[] storedLineagesAdded;
    protected List<Node>[] storedLineagesRemoved;
    
    //private List<Node>[] lineages;

    protected int intervalCount;
    protected int storedIntervalCount;
    
    private boolean lastIntervalDirty;

    /**
     * are the intervals known?
     */
    protected boolean intervalsKnown = false;
    protected double multifurcationLimit = -1.0;
    
    /**
     * =======================================================================
     * =======================================================================
     * =======================================================================
     * =======================================================================
     */
    
    public StructuredTreeIntervals() {
    }

    public StructuredTreeIntervals(Tree tree) {
//    	tree = 
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

    @SuppressWarnings("unchecked")
	@Override
    protected void restore() {
        intervalsKnown = false;
        Arrays.fill(intervalIsDirty, true);
        lastIntervalDirty = true;
        super.store();
        if (true)
        	return;

    	
        double[] tmp = storedIntervals;
        storedIntervals = intervals;
        intervals = tmp;

        int[] tmp2 = storedLineageCounts;
        storedLineageCounts = lineageCounts;
        lineageCounts = tmp2;

        int tmp3 = storedIntervalCount;
        storedIntervalCount = intervalCount;
        intervalCount = tmp3;
        
    	lineagesAdded = new List[intervalCount];
    	for (int i = 0; i < intervalCount; i++) {
    		if (storedLineagesAdded[i] != null) {
    			List<Node> nodeList = new ArrayList<>();
    			int nodeCount = storedLineagesAdded[i].size();
    			for (int n = 0; n < nodeCount; n++) {
    				nodeList.add(storedLineagesAdded[i].get(n).copy());
    			}
    			lineagesAdded[i] = nodeList;
    		}
    	}
    	
    	lineagesRemoved = new List[intervals.length];
    	for (int i = 0; i < intervalCount; i++) {
    		if (storedLineagesRemoved[i] != null) {
	    		List<Node> nodeList = new ArrayList<>();
	    		int nodeCount = storedLineagesRemoved[i].size();
	    		for (int n = 0; n < nodeCount; n++) {
	    			nodeList.add(storedLineagesRemoved[i].get(n).copy());
	    		}
	    		lineagesRemoved[i] = nodeList;
    		}
    	}


        
//       // Added these last two for struct coalescent models
//        lineagesAdded = storedLineagesAdded;        
//        lineagesRemoved = storedLineagesRemoved;
        
        setIntervalsClean();
        super.restore();
        
    }
    
    public void print(){
    	System.out.print("\n");
        System.out.println("restore tree intervals");
        for (int i = 0; i < 100; i++){
        	System.out.print("[");
        	if (lineagesRemoved[i]!=null)
        		for (int j = 0; j < lineagesRemoved[i].size(); j++)
        			System.out.print(lineagesRemoved[i].get(j).getNr() + " ,");
        	
        	System.out.print(lineagesAdded[i].get(0).getNr() + " ,");
        	System.out.print("]");
       }
    	System.out.print("\n");
    }

    @Override
    protected void store() {
        intervalsKnown = false;
        Arrays.fill(intervalIsDirty, true);
        lastIntervalDirty = true;
        super.restore();
        if (true)
        	return;

        // stores the lineage Counts per intervall and the intervalls in the arrays stored...
        System.arraycopy(lineageCounts, 0, storedLineageCounts, 0, lineageCounts.length);
        System.arraycopy(intervals, 0, storedIntervals, 0, intervals.length);
        storedIntervalCount = intervalCount;
        
        // Create new deep copies for storedLineagsAdded/Removed
        storedLineagesAdded = deepCopyLineagesAdded();
        storedLineagesRemoved = deepCopyLineagesRemoved();
       
        setIntervalsClean();
        super.store();
    }    
    
    @SuppressWarnings("unchecked")
	public List<Node>[] deepCopyLineagesAdded() {
    	List<Node>[] newList = new List[intervalCount];
    	for (int i = 0; i < intervalCount; i++) {
    		if (lineagesAdded[i] != null) {
    			List<Node> nodeList = new ArrayList<>();
    			int nodeCount = lineagesAdded[i].size();
    			for (int n = 0; n < nodeCount; n++) {
    				nodeList.add(lineagesAdded[i].get(n).copy());
    			}
    			newList[i] = nodeList;
    		}
    	}
    	return newList;
    	
    }
    
    @SuppressWarnings("unchecked")
	public List<Node>[] deepCopyLineagesRemoved() {    	
    	List<Node>[] newList = new List[intervals.length];
    	for (int i = 0; i < intervalCount; i++) {
    		if (lineagesRemoved[i] != null) {
	    		List<Node> nodeList = new ArrayList<>();
	    		int nodeCount = lineagesRemoved[i].size();
	    		for (int n = 0; n < nodeCount; n++) {
	    			nodeList.add(lineagesRemoved[i].get(n).copy());
	    		}
	    		newList[i] = nodeList;
	    		}
    	}
    	return newList;
    	
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
        return intervalCount;
    }

    /**
     * Gets an interval.
     */
    public double getInterval(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < 0 || i >= intervalCount){
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
    
    public double[] getCoalescentTimes(double[] coalescentTimes) {

        if (!intervalsKnown) {
            calculateIntervals();
        }

        if (coalescentTimes == null) coalescentTimes = new double[getSampleCount()];

        double time = 0;
        int coalescentIndex = 0;
        for (int i = 0; i < intervals.length; i++) {
            time += intervals[i];
            for (int j = 0; j < getCoalescentEvents(i); j++) {
                coalescentTimes[coalescentIndex] = time;
                coalescentIndex += 1;
            }
        }
        return coalescentTimes;
    }

    /**
     * Returns the number of uncoalesced lineages within this interval.
     * Required for s-coalescents, where new lineages are added as
     * earlier samples are come across.
     */
    public int getLineageCount(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        return lineageCounts[i];
    }
    
    public List<Node> getLineagesAdded(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        return lineagesAdded[i];
    }
    
    public List<Node> getLineagesRemoved(int index) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (index >= intervalCount) throw new IllegalArgumentException();
        return lineagesRemoved[index];
    }

    /**
     * Returns the number of coalescent events in an interval
     */
    public int getCoalescentEvents(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        if (i < intervalCount - 1) {
            return lineageCounts[i] - lineageCounts[i + 1];
        } else {
            return lineageCounts[i] - 1;
        }
    }

    /**
     * Returns the type of interval observed.
     */
    public IntervalType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        int numEvents = getCoalescentEvents(i);

        if (numEvents > 0) return IntervalType.COALESCENT;
        else if (numEvents < 0) return IntervalType.SAMPLE;
        else return IntervalType.NOTHING;
    }

    /**
     * @param i
     * @return intervalIsDirty returns if interval requires recalculation
     */
    protected boolean intervalIsDirty(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < intervalCount)
        	return intervalIsDirty[i];
        else
        	return lastIntervalDirty;
        	
    }    
    
    protected void setIntervalIsDirty(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < intervalCount){
        	intervalIsDirty[i] = true;
        }else{
        	lastIntervalDirty = true;
        }
    }    
  
    protected void setIntervalsClean(){
        intervalIsDirty = new boolean[treeInput.get().getNodeCount()];
        lastIntervalDirty = false;
    }


    /**
     * Recalculates all the intervals for the given beast.tree.
     */
    @SuppressWarnings("unchecked")
	public void calculateIntervals() {
//    	System.out.println("update intervals " + treeInput.get());
//    	if (tree!=null);
//    	else
		tree = treeInput.get();

		
//		System.out.println(String.format("%.200s",tree));
		
        final int nodeCount = tree.getNodeCount();

        double[] times = new double[nodeCount];
        int[] childCounts = new int[nodeCount];

        collectTimes(tree, times, childCounts);

        int[] indices = new int[nodeCount];        

        HeapSort.sort(times, indices);
        
        intervalIsDirty = new boolean[nodeCount];

        if (intervals == null || intervals.length != nodeCount) {
//        	System.out.println("dfllfdlk");
        	// Initialize new Lists
        	intervals = new double[nodeCount];
            lineageCounts = new int[nodeCount];
            lineagesAdded = new List[nodeCount];
            lineagesRemoved = new List[nodeCount];
            intervalIsDirty = new boolean[nodeCount];
            
            storedIntervals = new double[nodeCount];
            storedLineageCounts = new int[nodeCount];
            storedLineagesAdded = new List[nodeCount];
            storedLineagesRemoved = new List[nodeCount];
            
        } else {
            for (List<Node> l : lineagesAdded) {
            	if (l != null) l.clear();
            }
            for (List<Node> l : lineagesRemoved) {
            	if (l != null) l.clear();
            }
        }

        // start is the time of the first tip
        double start = times[indices[0]];
        int numLines = 0;
        int nodeNo = 0;
        intervalCount = 0;
        while (nodeNo < nodeCount) {

            int lineagesRemoved = 0;
            int lineagesAdded = 0;

            double finish = times[indices[nodeNo]];
            double next;

            do {
                final int childIndex = indices[nodeNo];
                final int childCount = childCounts[childIndex];
                
                
                // don't use nodeNo from here on in do loop
                nodeNo += 1;
                if (childCount == 0) {
                    addLineage(intervalCount, tree.getNode(childIndex));
                    lineagesAdded += 1;
                } else {
                    lineagesRemoved += (childCount - 1);

                    // record removed lineages
                    final Node parent = tree.getNode(childIndex);

                    // TODO avoid recalculation of one interval too many due
                    // to nodes that turned filthy when their parent changes
                    if (parent.isDirty()>0){
                    	intervalIsDirty[intervalCount] = true;
                    }else{
                    	intervalIsDirty[intervalCount] = false;                	
                    }
                    
                    
                    //assert childCounts[indices[nodeNo]] == beast.treeInput.get().getChildCount(parent);
                    //for (int j = 0; j < lineagesRemoved + 1; j++) {
                    for (int j = 0; j < childCount; j++) {
                        Node child = j == 0 ? parent.getLeft() : parent.getRight();
                        removeLineage(intervalCount, child);

                    }

                    // record added lineages
                    addLineage(intervalCount, parent);
                   
                    
                    // no mix of removed lineages when 0 th
                    if (multifurcationLimit == 0.0) {
                        break;
                    }
                }

                if (nodeNo < nodeCount) {
                    next = times[indices[nodeNo]];
                } else break;
                
            } while (Math.abs(next - finish) <= multifurcationLimit);

        	double diff = finish - start;
        	
        	if (diff != storedIntervals[intervalCount])
            	intervalIsDirty[intervalCount] = true;
       		
        	
        	if (lineagesAdded > 0) {

                if (intervalCount > 0 || ((finish - start) > multifurcationLimit)) {
                    intervals[intervalCount] = finish - start;
                    lineageCounts[intervalCount] = numLines;
                    intervalCount += 1;
                }

                start = finish;
            }

            // add sample event
            numLines += lineagesAdded;

            if (lineagesRemoved > 0) {        		

                intervals[intervalCount] = diff;
                lineageCounts[intervalCount] = numLines;
                intervalCount += 1;
                start = finish;
            }
            // coalescent event
            numLines -= lineagesRemoved;
        }    
//        print();
        intervalsKnown = true;
    }

    protected void addLineage(int interval, Node node) {
        if (lineagesAdded[interval] == null) lineagesAdded[interval] = new ArrayList<Node>();
        lineagesAdded[interval].add(node);
    }

    protected void removeLineage(int interval, Node node) {
        if (lineagesRemoved[interval] == null) lineagesRemoved[interval] = new ArrayList<Node>();
        lineagesRemoved[interval].add(node);
    }


    /**
     * extract coalescent times and tip information into array times from beast.tree.
     *
     * @param tree        the beast.tree
     * @param times       the times of the nodes in the beast.tree
     * @param childCounts the number of children of each node
     */
    protected static void collectTimes(Tree tree, double[] times, int[] childCounts) {
        Node[] nodes = tree.getNodesAsArray();
        for (int i = 0; i < nodes.length; i++) {
            Node node = nodes[i];
            times[i] = node.getHeight();
            childCounts[i] = node.isLeaf() ? 0 : 2;
        }
    }
          
    
    /**
     * get the total height of the genealogy represented by these
     * intervals.
     */
    public double getTotalDuration() {

        if (!intervalsKnown) {
            calculateIntervals();
        }
        double height = 0.0;
        for (int j = 0; j < intervalCount; j++) {
            height += intervals[j];
        }
        return height;
    }

    /**
     * Checks whether this set of coalescent intervals is fully resolved
     * (i.e. whether is has exactly one coalescent event in each
     * subsequent interval)
     */
    public boolean isBinaryCoalescent() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) > 0) {
                if (getCoalescentEvents(i) != 1) return false;
            }
        }

        return true;
    }

    /**
     * Checks whether this set of coalescent intervals coalescent only
     * (i.e. whether is has exactly one or more coalescent event in each
     * subsequent interval)
     */
    public boolean isCoalescentOnly() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) < 1) return false;
        }

        return true;
    }

    public double[] getIntervalStartEnd(int i){
    	double[] times = new double[2];
    	if (i == 0){
    		times[0] = 0.0;
    	}else{
    		List<Node> addedNodes = getLineagesAdded(i-1);
    		times[0] = addedNodes.get(0).getHeight();
    	}
		List<Node> addedNodes = getLineagesAdded(i);
		times[1] = addedNodes.get(0).getHeight();
		
		return times;   	
    }
    
    /**
     * Added Function
     */
    public void swap(){
    	calculateIntervals();
    	ArrayList<Integer> activeLineages = new ArrayList<Integer>();
        
        for (int i = 0; i < intervalCount; i++){
        	if(IntervalType.SAMPLE == getIntervalType(i)){
        		List<Node> incomingLines = getLineagesAdded(i);
        		for (Node l : incomingLines) {
        			activeLineages.add(l.getNr());
        		}
        	}

       	
        	if(IntervalType.COALESCENT == getIntervalType(i)){
        		List<Node> daughter = getLineagesRemoved(i);
            	
            	if (daughter.size() > 2) {
            		System.err.println("Multifurcation");
        		}
        		int d1 = activeLineages.indexOf(daughter.get(0).getNr());
        		int d2 = activeLineages.indexOf(daughter.get(1).getNr());
        		int j = i;
        		boolean swap = false;
        		
        		while (d1 == -1 || d2 == -1 && j < intervalCount){										// If true the nodes are in the wrong order
        			// Go to next event
        			j++;
                	if(IntervalType.COALESCENT == getIntervalType(j)){
                		daughter = getLineagesRemoved(j);
                		d1 = activeLineages.indexOf(daughter.get(0).getNr());
                		d2 = activeLineages.indexOf(daughter.get(1).getNr());
                	}
            		swap = true;
        		} 
        		
        		if(!swap){
	        		if (d1 > d2){
	        			activeLineages.remove(d1);
	        			activeLineages.remove(d2);
	        		}
	        		else{
	        			activeLineages.remove(d2);
	        			activeLineages.remove(d1);
	        		}
	        		List<Node> incomingLines = getLineagesAdded(i);
	        		activeLineages.add(incomingLines.get(0).getNr());
        		}
        		// Add parent
        		if (j == intervalCount){
        			break;
        		}       		
        		
        		if(swap){
        		    double inter = intervals[j];
        		    intervals[j] = intervals[i];
        		    intervals[i] = inter;        		    		
        		    double storedInter = storedIntervals[j];
        		    storedIntervals[j] = storedIntervals[i];
        		    storedIntervals[i] = storedInter;
        		    /**
        		     * The number of uncoalesced lineages within a particular interval.
        		     */
        		    /**
        		     * The lineages in each interval (stored by node ref).
        		     */
        		    List<Node> lineageAdded = lineagesAdded[j];
        		    lineagesAdded[j] = lineagesAdded[i];
        		    lineagesAdded[i] = lineageAdded;
        		    List<Node> lineageRemoved = lineagesRemoved[j];
        		    lineagesRemoved[j] = lineagesRemoved[i];
        		    lineagesRemoved[i] = lineageRemoved;
        			i--;
      			
        		}
        	}
        }
    }

}