package beast.mascot.empirical;

import org.jblas.util.Random;

import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class empiricalTreeDistributionOperator extends Operator {
    final public Input<empiricalTreeDistribution> empiricalTreeDistributionInput = new Input<>("empiricalTreeDistribution", "filename of tree file to be read", Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<Tree>("tree","");

    
    @Override
	public void initAndValidate() {
	}

	@Override
	public double proposal() {
		final Tree tree = treeInput.get(this);
		double treeNr =  Random.nextDouble()*empiricalTreeDistributionInput.get().trees.size();	
		Tree newTree = empiricalTreeDistributionInput.get().trees.get((int) treeNr);
		Node root = newTree.getRoot();
        root.setParent(null);
        tree.setRoot(root);

		
		return 0;
	}

}
