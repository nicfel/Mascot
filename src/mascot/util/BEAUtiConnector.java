package mascot.util;

import beast.base.core.BEASTInterface;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import mascot.dynamics.StructuredSkyline;
import mascot.parameterdynamics.NotSet;

/**
 * Class containing a static method used in the BEAUti template to tidy
 * up some loose ends.
 */
public class BEAUtiConnector {

    public static boolean customConnector(BeautiDoc doc) {
        for (BEASTInterface p : doc.getPartitions("tree")) {
            String pId = BeautiDoc.parsePartition(p.getID());
 
            StructuredSkyline dummy = (StructuredSkyline) doc.pluginmap.get("StructuredSkyline.t:" + pId);

            dummy.initAndValidate();
            
            if (dummy.parametricFunctionInput.get().neDynamicsInput.get().size()==dummy.getDimension()) {
            	// update all the names of the Ne dynamics are correct
                for (int i = 0; i < dummy.getDimension(); i++) {                	
                	if (!dummy.parametricFunctionInput.get().neDynamicsInput.get().get(i).getID().contentEquals(
                			"NeDynamics." + dummy.getStringStateValue(i) + ".t:"+pId)) {                		
                		dummy.parametricFunctionInput.get().neDynamicsInput.get().get(i).setID("NeDynamics." + dummy.getStringStateValue(i) + ".t:" + pId);
                		
                	}
                }
                
            	continue;
            }
            
            
            System.out.println("reset parameter");
            
            dummy.parametricFunctionInput.get().neDynamicsInput.get().clear();
            for (int i = 0; i < dummy.getDimension(); i++) {
            	NotSet notSet = new NotSet();
            	notSet.setID("NeDynamics." + dummy.getStringStateValue(i) + ".t:"+pId);
            	dummy.parametricFunctionInput.get().neDynamicsInput.get().add(notSet);
            }
            dummy.initAndValidate();
        }       

        return false;
    }
    
    public static boolean customConnectorPriorCleaner(BeautiDoc doc) {
    	CompoundDistribution prior = (CompoundDistribution) doc.pluginmap.get("prior");
        for (Distribution p : prior.pDistributions.get()) {
        	if (p instanceof Prior) {
	        	Prior pri = (Prior) p;
	        	
	        	if (pri.m_x.get() instanceof RealParameter) {
		        	RealParameter rp =  (RealParameter) pri.m_x.get();
	
	        		
		        	if (!rp.isEstimatedInput.get()) {
		        		System.out.println("----");
		        		System.out.println(p.getID());
		        		doc.disconnect(pri, "prior", "distribution");
		        	}
	        	}else if (pri.m_x.get() instanceof Difference) {
		        	RealParameter rp =  (RealParameter) ((Difference)  pri.m_x.get()).functionInput.get();
		        	if (!rp.isEstimatedInput.get()) {
		        		System.out.println("----");
		        		System.out.println(p.getID());
		        		doc.disconnect(pri, "prior", "distribution");
		        	}

	        	}else if (pri.m_x.get() instanceof First) {
		        	RealParameter rp =  (RealParameter) ((First)  pri.m_x.get()).functionInput.get();
		        	if (!rp.isEstimatedInput.get()) {
		        		System.out.println("----");
		        		System.out.println(p.getID());
		        		doc.disconnect(pri, "prior", "distribution");
		        	}

	        	}

        	}
        		
        }   


        return false;
    }
    
    

    
    
    
}
