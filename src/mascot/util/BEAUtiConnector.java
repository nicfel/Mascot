package mascot.util;

import beastfx.app.inputeditor.BeautiDoc;
import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.operator.UpDownOperator;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import mascot.dynamics.StructuredSkyline;
import mascot.parameterdynamics.NotSet;

import java.util.*;

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
            
            if (dummy.parametricFunctionInput.get().size()==dummy.getDimension()) {
            	//update all the names of the Ne dynamics are correct
                for (int i = 0; i < dummy.getDimension(); i++) {
                	String idName = dummy.parametricFunctionInput.get().get(i).getID().substring(
                			dummy.parametricFunctionInput.get().get(i).getID().indexOf(":") + 1,
                			dummy.parametricFunctionInput.get().get(i).getID().length());
                	if (!idName.contentEquals(dummy.getStringStateValue(i) + "."+pId)) {
                		String newID = dummy.parametricFunctionInput.get().get(i).getID().substring(0,
                				dummy.parametricFunctionInput.get().get(i).getID().indexOf(":") + 1);
                		dummy.parametricFunctionInput.get().get(i).setID(newID + "" + dummy.getStringStateValue(i) + "." + pId);
//                		 change all input id's
                		System.out.println(dummy.parametricFunctionInput.get().get(i).getInputs());
                		for (String s : dummy.parametricFunctionInput.get().get(i).getInputs().keySet()) {
                			System.out.println(dummy.parametricFunctionInput.get().get(i).getInputs().get(s).get());
                			if (dummy.parametricFunctionInput.get().get(i).getInputs().get(s).get()
                				instanceof Parameter){
                				BEASTObject obj = (BEASTObject) dummy.parametricFunctionInput.get().get(i).getInputs().get(s).get();
                				System.out.println("id is " + obj.getID());
                				
                				BEASTInterface bi = doc.pluginmap.get(obj.getID());

                				String old = bi.getID().substring(obj.getID().indexOf(":")+1,bi.getID().length());
                				System.out.println(bi.getID().replace(old, dummy.getStringStateValue(i) + "." + pId));
                				bi.setID(bi.getID().replace(old, dummy.getStringStateValue(i) + "." + pId));
                				System.out.println("new id is " + bi.getID());

            				}

                		}
                		
                	}
                }
                
            	continue;
            }
            
            
            System.out.println("reset parameter");
            
            dummy.parametricFunctionInput.get().clear();
            for (int i = 0; i < dummy.getDimension(); i++) {
            	NotSet notSet = new NotSet();
            	notSet.setID("notset.t:" +dummy.getStringStateValue(i) + "."+pId);
//            	notSet.setID("notset.t:" +pId);
            	dummy.parametricFunctionInput.get().add(notSet);
            }
            dummy.initAndValidate();
        }

        

        return false;
    }
}
