package beast.mascot.empirical;

import org.jblas.util.Random;

import beast.core.Input;
import beast.core.Operator;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class empiricalTreeDistributionOperator extends TreeOperator {
    final public Input<empiricalTreeDistribution> empiricalTreeDistributionInput = new Input<>("empiricalTreeDistribution", "filename of tree file to be read", Input.Validate.REQUIRED);
//    final public Input<Tree> treeInput = new Input<Tree>("tree","");

    
    @Override
	public void initAndValidate() {
	}

	@Override
	public double proposal() {
		final Tree tree = treeInput.get(this);
		double treeNr =  Random.nextDouble()*empiricalTreeDistributionInput.get().trees.size();	
//		System.out.println("\nakward tree operator at work " + (int) treeNr + "\n");
		Tree newTree = empiricalTreeDistributionInput.get().trees.get((int) treeNr);
        tree.startEditing(this);  // we change the tree

        for (int i = 0; i < tree.getNodeCount(); i++){
        	final Node node = tree.getNode(i);
        	node.setHeight(newTree.getNode(i).getHeight());
        	node.setID(newTree.getNode(i).getID());
        	if (!newTree.getNode(i).isLeaf()){
        		node.setLeft(newTree.getNode(i).getLeft());
        		node.setRight(newTree.getNode(i).getRight());
        	}
        	if (!newTree.getNode(i).isRoot()){
        		node.setParent(newTree.getNode(i).getParent());
        	}else{
        		node.setParent(null);
        	}
        	node.makeDirty(tree.IS_FILTHY);        	
        }        	

		
		return 0.0;
	}

}
