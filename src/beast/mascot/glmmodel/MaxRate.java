package beast.mascot.glmmodel;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.mascot.dynamics.GLM;
import beast.math.distributions.ParametricDistribution;

import org.apache.commons.math.MathException;


@Description("Produces prior (log) probability of value x." +
        "If x is multidimensional, the components of x are assumed to be independent, " +
        "so the sum of log probabilities of all elements of x is returned as the prior.")
public class MaxRate extends Distribution {
    final public Input<GLM> GLMStepwiseModelInput = new Input<>("GLMmodel", "glm model input");
    final public Input<ParametricDistribution> distInput = new Input<>("distr", "distribution used to calculate prior, e.g. normal, beta, gamma.", Validate.REQUIRED);
    final public Input<Boolean> migrationOnlyInput = new Input<>("migrationOnly", "put prior only on migration rates", false);
    final public Input<Boolean> NeOnlyInput = new Input<>("NeOnly", "put prior only on migration rates", false);

    /**
     * shadows distInput *
     */
    protected ParametricDistribution dist;

    
    @Override
    public void initAndValidate() {
        dist = distInput.get();
        calculateLogP();
    }

    @Override
    public double calculateLogP() {
    	Double[] mig = GLMStepwiseModelInput.get().getAllCoalescentRate();
		Double[] coal = GLMStepwiseModelInput.get().getAllBackwardsMigration();
    	
    	
    	RealParameter dCoal = new RealParameter(coal);
    	RealParameter dMig = new RealParameter(mig);
    	
    	logP = 0.0;
    	    	
    	if (migrationOnlyInput.get()){
	        logP += dist.calcLogP(dMig);
    	}else{
	        logP += dist.calcLogP(dCoal);
	        logP += dist.calcLogP(dMig);
    	}
        if (logP == Double.POSITIVE_INFINITY) {
            logP = Double.NEGATIVE_INFINITY;
        }
        return logP;
    }

    /**
     * return name of the parameter this prior is applied to *
     */
    public String getParameterName() {
//        if (m_x.get() instanceof BEASTObject) {
//            return ((BEASTObject) m_x.get()).getID();
//        }
        return "";
    }

    @Override
    public void sample(State state, Random random) {

        if (sampledFlag)
            return;

        sampledFlag = true;

        // Cause conditional parameters to be sampled
        sampleConditions(state, random);

        // sample distribution parameters
//        Function x = m_x.get();
//
//        Double[] newx;
//        try {
//            newx = dist.sample(1)[0];
//
//            if (x instanceof RealParameter) {
//                for (int i = 0; i < newx.length; i++) {
//                    ((RealParameter) x).setValue(i, newx[i]);
//                }
//            } else if (x instanceof IntegerParameter) {
//                for (int i = 0; i < newx.length; i++) {
//                    ((IntegerParameter) x).setValue(i, (int)Math.round(newx[i]));
//                }
//            }
//
//        } catch (MathException e) {
//            e.printStackTrace();
//            throw new RuntimeException("Failed to sample!");
//        }
    }

    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(dist.getID());
        return conditions;
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();

        String id = null;
//        if (m_x.get() != null && m_x.get() instanceof BEASTInterface) {
//            arguments.add(((BEASTInterface)m_x.get()).getID());
//        }

        return arguments;
    }
}
