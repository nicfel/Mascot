package beast.mascot.empirical;

import java.util.List;

import org.jblas.util.Random;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.Tree;

public class InitFromEmpiricalTreeDistribution extends Tree implements StateNodeInitialiser {
    final public Input<empiricalTreeDistribution> empiricalTreeDistributionInput = 
    		new Input<>("empiricalTreeDistribution", 
    				"filename of tree file to be read", Input.Validate.REQUIRED);

		
    @Override
    public void initStateNodes() {
		double treeNr =  Random.nextDouble()*empiricalTreeDistributionInput.get().trees.size();	
		Tree newTree = empiricalTreeDistributionInput.get().trees.get((int) treeNr);
		
//		System.out.println(newTree);

        if (m_initial.get() != null) {
            m_initial.get().assignFrom(newTree);
        }
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
        if (m_initial.get() != null) {
            stateNodes.add(m_initial.get());
        }
    }

}
