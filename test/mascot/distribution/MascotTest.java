package mascot.distribution;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.coalescent.Coalescent;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.inference.parameter.RealParameter;
import junit.framework.Assert;
import mascot.dynamics.Constant;
import org.junit.Test;


public class MascotTest  {

	@Test
	public void testMascotStructured(){
		
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
		alignment.initByName("sequence", sequence1,"sequence", sequence2,"sequence", sequence3,"sequence", sequence4, "sequence", sequence5);
		
		TaxonSet taxa = new TaxonSet();
		taxa.initByName("alignment", alignment);
		
		//build trait set
		TraitSet traitSet = new TraitSet();
		traitSet.initByName("value", "a1=a,a2=a,a3=a,b1=b,b2=b", "traitname", "type", "taxa", taxa);
		
		Tree tree = new TreeParser("(((a1:1,a2:2):1,(b1:1,b2:1.5):2):1,a3:1)");
		Tree traitTree = new Tree();
		tree.initByName("taxonset", taxa, "trait", traitSet);
		StructuredTreeIntervals st = new StructuredTreeIntervals();
		st.initByName("tree",tree);
		
		// build constant dynamics
		RealParameter Ne = new RealParameter("1 2");
		RealParameter backwardsMigration = new RealParameter("0.3 2");
		int dim = 2;
		Constant constant = new Constant();
		constant.initByName("backwardsMigration", backwardsMigration, "Ne", Ne, "dimension", dim, "typeTrait", traitSet);

		
		Mascot mascot = new Mascot();
		
		mascot.initByName("structuredTreeIntervals", st, "dynamics", constant, "tree", tree);
		
		double logP = mascot.calculateLogP();
		Assert.assertEquals(logP,-6.870390751933608, 1e-15);
		

	}
	
	@Test
	public void testMascotUnstructured(){
		
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
		alignment.initByName("sequence", sequence1,"sequence", sequence2,"sequence", sequence3,"sequence", sequence4, "sequence", sequence5);
		
		
		TaxonSet taxa = new TaxonSet();
		taxa.initByName("alignment", alignment);
		
		//build trait set
		TraitSet traitSet = new TraitSet();
		traitSet.initByName("value", "a1=a,a2=a,a3=a,b1=a,b2=a", "traitname", "type", "taxa", taxa);
		
		Tree tree = new TreeParser("(((a1:1,a2:2):1,(b1:1,b2:1.5):2):1,a3:1)");
		Tree traitTree = new Tree();
		tree.initByName("taxonset", taxa, "trait", traitSet);
		StructuredTreeIntervals st = new StructuredTreeIntervals();
		st.initByName("tree",tree);
		
		// build constant dynamics
		RealParameter Ne = new RealParameter("1.5 2");
		RealParameter backwardsMigration = new RealParameter("0 0");
		int dim = 2;
		Constant constant = new Constant();
		constant.initByName("backwardsMigration", backwardsMigration, "Ne", Ne, "dimension", dim, "typeTrait", traitSet);

		
		Mascot mascot = new Mascot();
		
		mascot.initByName("structuredTreeIntervals", st, "dynamics", constant, "tree", tree);
		
		// also set up the constant coalescent
		TreeIntervals ti = new TreeIntervals();
		ti.initByName("tree",tree);
		
		RealParameter Nec = new RealParameter("1.5");
		
		
		ConstantPopulation cp = new ConstantPopulation();
		cp.initByName("popSize", Nec);
		
		Coalescent coalescent = new Coalescent();
		coalescent.initByName("populationModel", cp, "treeIntervals", ti);
		
		Assert.assertTrue(Math.abs(mascot.calculateLogP()-coalescent.calculateLogP())<0.000000001);
	}

}
