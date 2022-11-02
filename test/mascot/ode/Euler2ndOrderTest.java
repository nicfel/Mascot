package mascot.ode;

import java.util.Arrays;

import org.junit.Test;

import junit.framework.Assert;
import mascot.ode.Euler2ndOrder;


public class Euler2ndOrderTest  {
	
	@Test
	public void testValues(){
		int states = 2;
		
		// migration rates
		double[] migration_rates = new double[2*2];
		migration_rates[0*2+0] = 0;
		migration_rates[1*2+1] = 0;
		migration_rates[0*2+1] = 0.2;
		migration_rates[1*2+0] = 2;
		
		// coal rates
		double[] coalescent_rates = {1,2};
		
		// lins
		int lineages = 2;
		
		// lins
		double[] p = {0.5, 0.5, 1, 0, 1, 0, 0};				
		
		
		Euler2ndOrder euler2ndOrder =  new Euler2ndOrder(migration_rates, coalescent_rates, lineages+1, states, 0.001, 0.2);
				
		double[] pDot = new double[p.length];
		double[] pDotDot = new double[p.length];
		double[] pDotDotDot = new double[p.length];
		euler2ndOrder.calculateValues(1.0, p, pDot, pDotDot, pDotDotDot, p.length);
		
		// check that it does not depend if multiplicators are used
		Assert.assertTrue(Math.abs(p[0]-0.7280934851118206)<0.0000000001);
		Assert.assertTrue(Math.abs(p[1]-0.2719065148881794)<0.0000000001);
		Assert.assertTrue(Math.abs(p[2]-0.8961702265249307)<0.0000000001);
		Assert.assertTrue(Math.abs(p[3]-0.10382977347506926)<0.0000000001);
		Assert.assertTrue(Math.abs(p[4]-0.8961702265249307)<0.0000000001);		
		Assert.assertTrue(Math.abs(p[5]-0.10382977347506926)<0.0000000001);
		Assert.assertTrue(Math.abs(p[6]+2.1375326433098785)<0.0000000001);
	}

	@Test
	public void testDifferentSpecifications(){
		int states = 2;
		
		// mult
		int[] multiplicator = new int[2];
		multiplicator[0] = 1;
		multiplicator[1] = 2;
		
		// migration rates
		double[] migration_rates = new double[2*2];
		migration_rates[0*2+0] = 0;
		migration_rates[1*2+1] = 0;
		migration_rates[0*2+1] = 0.2;
		migration_rates[1*2+0] = 2;
		
		// coal rates
		double[] coalescent_rates = {1,2};
		
		// lins
		int lineages = 2;
		
		// lins
		double[] p = {0.5, 0.5, 1, 0, 1, 0, 0};		
		double[] p_mul = {0.5, 0.5, 1, 0, 0};
		
		
		
		Euler2ndOrder euler2ndOrder =  new Euler2ndOrder(migration_rates, coalescent_rates, lineages+1, states, 0.001, 0.2);
		Euler2ndOrder euler2ndOrderMult =  new Euler2ndOrder(multiplicator, migration_rates, coalescent_rates, lineages, states, 0.001, 0.2);
				
		double[] pDot = new double[p.length];
		double[] pDotDot = new double[p.length];
		double[] pDotDotDot = new double[p.length];
		euler2ndOrder.calculateValues(1.0, p, pDot, pDotDot, pDotDotDot, p.length);
		
		pDot = new double[p_mul.length];
		pDotDot = new double[p_mul.length];
		pDotDotDot = new double[p_mul.length];
		euler2ndOrderMult.calculateValues(1.0, p_mul, pDot, pDotDot, pDotDotDot, p_mul.length);
		
		// check that it does not depend if multiplicators are used
		Assert.assertTrue(Math.abs(p[0]-p_mul[0])<0.0000000001);
		Assert.assertTrue(Math.abs(p[4]-p_mul[2])<0.0000000001);		
		Assert.assertFalse(Math.abs(p[1]-p_mul[3])<0.0000000001);
		Assert.assertTrue(Math.abs(p[6]-p_mul[4])<0.0000000001);
		
		// check that it does not depend if indicators
		
		int[] indicators = new int[1*2];
		indicators[0*2+0] = 0;
		indicators[0*2+1] = 1;
		
		double[] p_ind = {0.5, 0.5, 1, 0, 1, 0, 0};
		double[] p_new = {0.5, 0.5, 1, 0, 1, 0, 0};


		pDot = new double[p.length];
		pDotDot = new double[p.length];
		pDotDotDot = new double[p.length];
		
		
		Euler2ndOrder euler2ndOrderInd = new Euler2ndOrder(migration_rates, indicators, coalescent_rates, lineages+1, states, 0.001, 0.2);
		
		migration_rates[0*2+0] = 0;
		migration_rates[1*2+1] = 0;
		migration_rates[0*2+1] = 0.2;
		migration_rates[1*2+0] = 0;
		
		euler2ndOrderInd.calculateValues(1.0, p_ind, pDot, pDotDot, pDotDotDot, p_ind.length);
		

		pDot = new double[p.length];
		pDotDot = new double[p.length];
		pDotDotDot = new double[p.length];		
		
		
		euler2ndOrder = new Euler2ndOrder(migration_rates, coalescent_rates, lineages+1, states, 0.001, 0.2);
		euler2ndOrder.calculateValues(1.0, p_new, pDot, pDotDot, pDotDotDot, p_new.length);
				
		Assert.assertTrue(Math.abs(p_new[0]-p_ind[0])<0.0000000001);
		Assert.assertTrue(Math.abs(p_new[1]-p_ind[1])<0.0000000001);		
		Assert.assertFalse(Math.abs(p_new[1]-p_ind[3])<0.0000000001);
		Assert.assertTrue(Math.abs(p_new[6]-p_ind[6])<0.0000000001);


	}

}
