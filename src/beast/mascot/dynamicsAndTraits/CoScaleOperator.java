/*
* File ScaleOperator.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
package beast.mascot.dynamicsAndTraits;


import java.text.DecimalFormat;
import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;


@Description("Scales a parameter or a complete beast.tree (depending on which of the two is specified.")
public class CoScaleOperator extends Operator {

    public final Input<RealParameter> parameterInput = new Input<>("parameter", "if specified, this parameter is scaled",
            Input.Validate.REQUIRED);

    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: larger means more bold proposals", 1.0);

    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 1.0 - 1e-8);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 1e-8);
    
    public Input<variableTrait> variableTraitInput = new Input<>("variableTrait",
            "traits that can change");


    /**
     * shadows input *
     */
    private double m_fScaleFactor;

    private double upper, lower;

    @Override
    public void initAndValidate() {
        m_fScaleFactor = scaleFactorInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();
    }


    protected boolean outsideBounds(final double value, final RealParameter param) {
        final Double l = param.getLower();
        final Double h = param.getUpper();

        return (value < l || value > h);
        //return (l != null && value < l || h != null && value > h);
    }

    protected double getScaler() {
        return (m_fScaleFactor + (Randomizer.nextDouble() * ((1.0 / m_fScaleFactor) - m_fScaleFactor)));
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	// get the trait mapping to know which indices to change at the same time
    	int[] minStateMapping = variableTraitInput.get().getMinimalStateMapping(true);
    	int parameterDimension = variableTraitInput.get().getNrStates();

    	
        try {

            double hastingsRatio;
            final double scale = getScaler();

            final RealParameter param = parameterInput.get(this);

            assert param.getLower() != null && param.getUpper() != null;

//            final int dim = param.getDimension();

            hastingsRatio = -Math.log(scale);

            // which position to scale
            final int index;

            // any is good
            index = Randomizer.nextInt(parameterDimension);            
            
            double oldValue = 0.0;
            for (int i = 0; i < minStateMapping.length; i++){
            	if (minStateMapping[i] == index){
            		oldValue = param.getValue(i);
            	}
            }

            if (oldValue == 0) {
                // Error: parameter has value 0 and cannot be scaled
                return Double.NEGATIVE_INFINITY;
            }

            final double newValue = scale * oldValue;

            if (outsideBounds(newValue, param)) {
                // reject out of bounds scales
                return Double.NEGATIVE_INFINITY;
            }
            
            // update parameter
            for (int i = 0; i < minStateMapping.length; i++){
            	if (minStateMapping[i] == index){
            		param.setValue(i, newValue);
            	}
            }

            return hastingsRatio;

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }


    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(1.0 / m_fScaleFactor - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return m_fScaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        m_fScaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(m_fScaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

} // class ScaleOperator
