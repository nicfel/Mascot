package beast.mascot.empirical;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jblas.util.Random;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

public class initFromTreeNr extends TreeParser implements StateNodeInitialiser {
    final public Input<String> treeFilenameInput = new Input<>("filename", "filename of tree file to be read", Input.Validate.REQUIRED);
    final public Input<Integer> treeNumberInput = new Input<>("treeNumber", "filename of tree file to be read", 1);

    private Tree tree;
    
    @Override
    public void initAndValidate() {
    	// read in the file
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(treeFilenameInput.get()));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Input file not found.");
        }
        String line;
        StringBuilder newickBuilder = new StringBuilder();
        
        // get translation
        ArrayList<String> translation1 = new ArrayList<String>();
        ArrayList<String> translation2 = new ArrayList<String>();
              
        boolean stop = false;
        int counter = 1;
        
        System.out.print("read in tree number " + treeNumberInput.get() +  " from file named " +  treeFilenameInput.get() + "\n");
        try {
            while ((line = reader.readLine()) != null && !stop){
            	// check for the first line of translate
            	if (line.trim().equals("Translate")){
            		do{
            			line = reader.readLine();
            			String[] tmp_line = line.trim().split(" "); 
            			if (tmp_line.length==2){
            				translation1.add(tmp_line[0]);
            				translation2.add(tmp_line[1].replaceAll(",", ""));
            			}
            		}while (!line.trim().equals(";"));        		
            	}
    			String[] tmp_line = line.trim().split(" "); 
    			
    			
    			while (tmp_line[0].equals("tree") && !stop){
    				if (counter == (int) treeNumberInput.get()){
	    				tree = getTree(tmp_line[tmp_line.length-1], translation1, translation2).copy();
	    				stop = true;
    				}
    				counter++;
        			line = reader.readLine();
        			if(line == null)
        				break;
        			tmp_line = line.trim().split(" ");
        			if (tmp_line.length<2)
        				break;        			
    			}
            }
        } catch (IOException e) {
            throw new RuntimeException("Error reading from input file.");
        }
       
        System.out.print("\nfinished reading tree");
    }
    
    private Tree getTree(String tree_string, ArrayList<String> translation1,  ArrayList<String> translation2){
    	//
    	for (int i = 0; i < translation1.size(); i++){
    		tree_string = tree_string.replace(String.format("(%s:", translation1.get(i)),
    				String.format("(%s:", translation2.get(i)));
    		tree_string = tree_string.replace(String.format(",%s:", translation1.get(i)),
    				String.format(",%s:", translation2.get(i)));
    	}
    	newickInput.setValue(tree_string, this);
        isLabelledNewickInput.setValue(true, this);
        adjustTipHeightsInput.setValue(false, this);
        allowSingleChildInput.setValue(false, this);
        offsetInput.setValue(1, this);
    	super.initAndValidate();
    	
    	Node root = parseNewick(tree_string);
    	return new Tree(root);
    }


    
		
//    @Override
//    public void initStateNodes() {
//    	System.out.println(tree);
//        if (m_initial.get() != null) {
//            m_initial.get().assignFrom(tree);
//        }
//    }
//
//    @Override
//    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
//        if (m_initial.get() != null) {
//            stateNodes.add(m_initial.get());
//        }
//    }

}
