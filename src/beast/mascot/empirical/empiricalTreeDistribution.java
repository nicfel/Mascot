package beast.mascot.empirical;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;


@Description("Read tree distribution from file, similar to the empiricalTreeDistributionOperator from BEAST1" +
			"and from TreeFromNewickFile from Tim Vaughan")
public class empiricalTreeDistribution extends TreeParser implements StateNodeInitialiser {
    final public Input<String> filenameInput = new Input<>("filename", "filename of tree file to be read", "tree.trees");

    ArrayList<Tree> trees;
    
    @Override
    public void initAndValidate() {
    	// read in the file
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(filenameInput.get()));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Input file not found.");
        }
        String line;
        StringBuilder newickBuilder = new StringBuilder();
        
        // get translation
        ArrayList<String> translation1 = new ArrayList<String>();
        ArrayList<String> translation2 = new ArrayList<String>();
        
        
        trees = new ArrayList<>();
        System.out.print("read in trees from file\n");
        try {
            while ((line = reader.readLine()) != null){
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
    			
    			
    			while (tmp_line[0].equals("tree")){
    		        System.out.print(".");
    				trees.add(getTree(tmp_line[tmp_line.length-1], translation1, translation2));
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
       
        System.out.print("\nfinished reading " + trees.size() + " trees");
    }
    
    private Tree getTree(String tree_string, ArrayList<String> translation1,  ArrayList<String> translation2){
    	//
    	for (int i = 0; i < translation1.size(); i++){
    		tree_string = tree_string.replace(String.format("(%s:", translation1.get(i)),
    				String.format("(%s:", translation2.get(i)));
    		tree_string = tree_string.replace(String.format(",%s:", translation1.get(i)),
    				String.format(",%s:", translation2.get(i)));
    	}
//    	System.out.println(tree_string);
//    	System.exit(0);
    	newickInput.setValue(tree_string, this);
        isLabelledNewickInput.setValue(true, this);
        adjustTipHeightsInput.setValue(false, this);
        allowSingleChildInput.setValue(false, this);
        offsetInput.setValue(1, this);
    	super.initAndValidate();
    	
    	Node root = parseNewick(tree_string);
//    	Tree tree = new Tree(root);
//    	System.out.println(tree);
    	
    	return new Tree(root);
    }

    /*
     *StateNodeInitializer implementation
     */
    @Override
    public void initStateNodes() { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodeList) {
        stateNodeList.add(this);
    }    


}
