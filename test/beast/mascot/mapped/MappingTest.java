package beast.mascot.mapped;

import java.util.Arrays;

import org.junit.Test;

import beast.app.mascot.beauti.TreeWithTrait;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.Coalescent;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.mascot.distribution.MappedMascot;
import beast.mascot.distribution.Mascot;
import beast.mascot.distribution.StructuredTreeIntervals;
import beast.mascot.dynamics.Constant;
import beast.mascot.logger.StructuredTreeLogger;
import beast.mascot.ode.Euler2ndOrder;
import beast.util.Randomizer;
import beast.util.TreeParser;
import coalre.network.Network;
import junit.framework.Assert;
import score.distribution.SCORE;
import score.dynamics.ConstantReassortment;


public class MappingTest  {

	
	@Test
	public void testMapping2() throws Exception{
		
		//build alignment
		Sequence sequence1 = new Sequence();
		Sequence sequence2 = new Sequence();
		Sequence sequence3 = new Sequence();
		Sequence sequence4 = new Sequence();
		Sequence sequence5 = new Sequence();
		sequence1.initByName("taxon", "a1", "value", "???");
		sequence2.initByName("taxon", "a2", "value", "???");
		sequence3.initByName("taxon", "a3", "value", "???");
		sequence4.initByName("taxon", "b1", "value", "???");
		sequence5.initByName("taxon", "b2", "value", "???");
		
		
		Alignment alignment = new Alignment();
		alignment.initByName("sequence", sequence1,"sequence", sequence2,"sequence", 
				sequence3,"sequence", sequence4,"sequence", sequence5);
		
		TaxonSet taxa = new TaxonSet();
		taxa.initByName("alignment", alignment);
		
		//build trait set
		TraitSet traitSet = new TraitSet();
		traitSet.initByName("value", "a1=a,a2=a,a3=a,b1=b,b2=b", "traitname", "type", "taxa", taxa);
		
		Tree tree = new TreeParser("(((a1:0.5,a2:0.1):0.5,(b2:0.1,b1:0.1):0.2):.01,a3:0.01)");
		
		System.out.println(tree);
		
		Tree traitTree = new Tree();
		tree.initByName("taxonset", taxa, "trait", traitSet);
		StructuredTreeIntervals st = new StructuredTreeIntervals();
		st.initByName("tree",tree);
		
		// build constant dynamics
		RealParameter Ne = new RealParameter("1 0.1");
		RealParameter backwardsMigration = new RealParameter("0.5 1");
		int dim = 2;
		Constant constant = new Constant();
		constant.initByName("backwardsMigration", backwardsMigration, "Ne", Ne, "dimension", dim, "typeTrait", traitSet);

		
		Mascot mascot = new Mascot();
		
		mascot.initByName("structuredTreeIntervals", st, "dynamics", constant, "tree", tree);
		
		// compute the node probs
		StructuredTreeLogger stl = new StructuredTreeLogger();
		stl.initByName("mascot", mascot, "epsilon", 0.0000001, "upDown", true);
		
		StructuredTreeLogger stl2 = new StructuredTreeLogger();
		stl2.initByName("mascot", mascot, "epsilon", 0.0000001, "upDown", false);

		
		MappedMascot mappedMascot = new MappedMascot();
		mappedMascot.initByName("structuredTreeIntervals", st, 
				"dynamics", constant, 
				"tree", tree,
				"maxIntegrationStepMapping",0.001);
	
		
		stl.calcForTest();
		System.out.println();
		stl2.calcForTest();
		
		Randomizer.setSeed(4);
		
		
		Node[] nodes = tree.getNodesAsArray();
		
		double[][] counts = new double[nodes.length][dim];
		int nrreps = 10000;
		
		
		for (int r = 0; r < nrreps; r++) {
//			System.out.println();
			mappedMascot.calculateLogP();
			for (int i = 0; i < nodes.length;i++) {
				if (!nodes[i].isLeaf()) {
					counts[nodes[i].getNr()][mappedMascot.getNodeState(nodes[i].getNr())]++;
				}
			}
		}
 
		for (int a = 0; a < counts.length; a++)
			for (int b = 0; b < counts[a].length; b++)
				counts[a][b] /= nrreps;
		

		
		double diff = 0.0;
		
		for (int i = 0; i < nodes.length;i++) {
			if (!nodes[i].isLeaf()) {
				double[] nodeProbs = stl.getStateProbOnly(nodes[i].getNr()).toArray();
				double[] nodeProbs2 = stl2.getStateProbOnly(nodes[i].getNr()).toArray();
				diff += Math.abs(counts[nodes[i].getNr()][0] - nodeProbs[0]);
				
				System.out.println(Arrays.toString(counts[nodes[i].getNr()]));
				System.out.println(Arrays.toString(nodeProbs));
				System.out.println(Arrays.toString(nodeProbs2));
				System.out.println(nodes[i].getNr());
				System.out.println(nodes[i].getHeight());
				System.out.println(Math.abs(counts[nodes[i].getNr()][0] - nodeProbs[0]));
				System.out.println("");
			}
		}
		
		
		mappedMascot.log(0, System.out);
		System.out.println();
		stl.log(0, System.out);
		System.out.println();
		
		System.out.println(tree);
		System.exit(0);

	}

	
	
	@Test
	public void testMapping() throws Exception{
		
		//build alignment
		Sequence sequence1 = new Sequence();
		Sequence sequence2 = new Sequence();
		Sequence sequence3 = new Sequence();
		Sequence sequence4 = new Sequence();
		Sequence sequence5 = new Sequence();
		Sequence sequence6 = new Sequence();
		Sequence sequence7 = new Sequence();
		sequence1.initByName("taxon", "a1", "value", "???");
		sequence2.initByName("taxon", "a2", "value", "???");
		sequence3.initByName("taxon", "a3", "value", "???");
		sequence4.initByName("taxon", "b1", "value", "???");
		sequence5.initByName("taxon", "b2", "value", "???");
		sequence6.initByName("taxon", "c1", "value", "???");
		sequence7.initByName("taxon", "c2", "value", "???");
		
		
		Alignment alignment = new Alignment();
		alignment.initByName("sequence", sequence1,"sequence", sequence2,"sequence", 
				sequence3,"sequence", sequence4, "sequence", sequence5,
				"sequence",sequence6, "sequence", sequence7);
		
		TaxonSet taxa = new TaxonSet();
		taxa.initByName("alignment", alignment);
		
		//build trait set
		TraitSet traitSet = new TraitSet();
		traitSet.initByName("value", "a1=a,a2=a,a3=a,b1=b,b2=b,c1=a,c2=b", "traitname", "type", "taxa", taxa);
		
		Tree tree = new TreeParser("((((a1:0.5,a2:0.1):0.5,(b1:0.5,b2:1):0.2):1,a3:1):0.1,(c1:0.5,c2:1.5):0.3)");
		Tree traitTree = new Tree();
		tree.initByName("taxonset", taxa, "trait", traitSet);
		StructuredTreeIntervals st = new StructuredTreeIntervals();
		st.initByName("tree",tree);
		
		// build constant dynamics
		RealParameter Ne = new RealParameter("1 1");
		RealParameter backwardsMigration = new RealParameter("1 1");
		int dim = 2;
		Constant constant = new Constant();
		constant.initByName("backwardsMigration", backwardsMigration, "Ne", Ne, "dimension", dim, "typeTrait", traitSet);

		
		Mascot mascot = new Mascot();
		
		mascot.initByName("structuredTreeIntervals", st, "dynamics", constant, "tree", tree);
		
		// compute the node probs
		StructuredTreeLogger stl = new StructuredTreeLogger();
		stl.initByName("mascot", mascot, "epsilon", 0.0000001, "upDown", true);
		
		
		MappedMascot mappedMascot = new MappedMascot();
		mappedMascot.initByName("structuredTreeIntervals", st, 
				"dynamics", constant, 
				"tree", tree,
				"maxIntegrationStepMapping",0.001);
	
		
		stl.calcForTest();
		
		Randomizer.setSeed(4);
		
		
		Node[] nodes = tree.getNodesAsArray();
		
		double[][] counts = new double[nodes.length][dim];
		int nrreps = 100000;
		
		
		for (int r = 0; r < nrreps; r++) {
//			System.out.println();
			mappedMascot.calculateLogP();
			for (int i = 0; i < nodes.length;i++) {
				if (!nodes[i].isLeaf()) {
					counts[nodes[i].getNr()][mappedMascot.getNodeState(nodes[i].getNr())]++;
				}
			}
		}
//		System.exit(0);
 
		for (int a = 0; a < counts.length; a++)
			for (int b = 0; b < counts[a].length; b++)
				counts[a][b] /= nrreps;
		

		
		double diff = 0.0;
		
		for (int i = 0; i < nodes.length;i++) {
			if (!nodes[i].isLeaf()) {
				double[] nodeProbs = stl.getStateProbOnly(nodes[i].getNr()).toArray();
				diff += Math.abs(counts[nodes[i].getNr()][0] - nodeProbs[0]);
				
				System.out.println(Arrays.toString(counts[nodes[i].getNr()]));
				System.out.println(Arrays.toString(nodeProbs));
				System.out.println(nodes[i].getNr());
				System.out.println(nodes[i].getHeight());
				System.out.println(Math.abs(counts[nodes[i].getNr()][0] - nodeProbs[0]));
				System.out.println("");
			}
		}
		
		
//		System.out.println();
//		System.out.println(diff/6);

		
		mappedMascot.log(0, System.out);
		System.out.println();
		stl.log(0, System.out);
		System.out.println();
		
		System.out.println(tree);

//		traitSet.initByName("value", "a1=a,a2=a,a3=a,b1=b,b2=b,c1=c,c2=c", "traitname", "type", "taxa", taxa);
//		
//		Tree tree = new TreeParser("((((a1:1,a2:2):1,(b1:1,b2:1.5):2):1,a3:1):0.1,(c1:0.5,c2:1.5):0.3)");
//		Tree traitTree = new Tree();
//		tree.initByName("taxonset", taxa, "trait", traitSet);
//		StructuredTreeIntervals st = new StructuredTreeIntervals();
//		st.initByName("tree",tree);
//		
//		// build constant dynamics
//		RealParameter Ne = new RealParameter("0.5 2 1");
//		RealParameter backwardsMigration = new RealParameter("0.1 0.3 0.2 0.1 0.1 0.4");
//		int dim = 3;
//		Constant constant = new Constant();
//		constant.initByName("backwardsMigration", backwardsMigration, "Ne", Ne, "dimension", dim, "typeTrait", traitSet);
		
		// sample from mapping

	}
	
	



}
