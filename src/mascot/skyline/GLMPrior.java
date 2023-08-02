package mascot.skyline;

import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import mascot.glmmodel.CovariateList;
import mascot.parameterdynamics.NeDynamics;
import mascot.parameterdynamics.NeDynamicsList;
import mascot.parameterdynamics.Skygrowth;

import java.util.List;
import java.util.Random;

public class GLMPrior extends Distribution {

    public Input<CovariateList> covariateListInput = new Input<>("covariateList", "input of covariates", Input.Validate.REQUIRED);
    public Input<RealParameter> scalerInput = new Input<>("scaler", "input of covariates scaler", Input.Validate.REQUIRED);
    public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "input of covariates scaler", Input.Validate.REQUIRED);

    public Input<RateShifts> rateShiftsInput = new Input<>(
            "rateShifts", "input of timings of rate shifts relative to the most recent sample", Input.Validate.OPTIONAL);
    public Input<NeDynamicsList> NeFunctionInput = new Input<>(
            "NeDynamics", "input of the log effective population sizes", Input.Validate.REQUIRED);

    final public Input<ParametricDistribution> distInput = new Input<>("distr",
            "distribution used to on the error terms of the GLM.",
            Input.Validate.REQUIRED);

    protected ParametricDistribution dist;

    double[] intTimes;

    int firstlargerzero;

    double[] differences;
    double mean;

    int verticalEntries;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        rateShiftsInput.get().initAndValidate();

        // set the dimension of the scalers, indicators and potentially the error term
        scalerInput.get().setDimension(covariateListInput.get().size());
        indicatorInput.get().setDimension(covariateListInput.get().size());



        // set the number of intervals for the GLM models
//        if (fromBeautiInput.get()) {
//            // the type to trait map is needed to read in predictors
//            migrationGLMInput.get().covariateListInput.get().traitToType = new HashMap<>(traitToType);
//            NeGLMInput.get().covariateListInput.get().traitToType = new HashMap<>(traitToType);
//
//            migrationGLMInput.get().covariateListInput.get().nrIntervals = rateShiftsInput.get().getDimension();
//            NeGLMInput.get().covariateListInput.get().nrIntervals = rateShiftsInput.get().getDimension();
//
//            migrationGLMInput.get().setNrDummy();
//            NeGLMInput.get().setNrDummy();
//        }
//        if (rateShiftsInput.get().getDimension()>0) {
//            if (GLMInput.get().covariateListInput.get().size()>0)
//                GLMInput.get().setNrIntervals(rateShiftsInput.get().getDimension(), dimensionInput.get(), true);
//        }

        //get the first non zero element
        firstlargerzero = 0;
        for (int i=0; i < rateShiftsInput.get().getDimension(); i++) {
            if (rateShiftsInput.get().getValue(i)>0) {
                firstlargerzero=i;
                break;
            }
        }
        // initialize the intervals
        intTimes = new double[rateShiftsInput.get().getDimension()-firstlargerzero];
        for (int i=0; i < intTimes.length; i++) {
            if (i==0) {
                intTimes[i] = rateShiftsInput.get().getValue(i+firstlargerzero);
            }
            else {
                intTimes[i] = rateShiftsInput.get().getValue(i+firstlargerzero)-rateShiftsInput.get().getValue(i-1+firstlargerzero);
            }
        }

        verticalEntries = covariateListInput.get().get(0).getDimension()/(rateShiftsInput.get().getDimension());

    }


    public double calculateLogP() {
        logP = 0;

        // compute all the error terms
        NeDynamics ne = NeFunctionInput.get().neDynamicsInput.get().get(0);
        if (ne instanceof Skygrowth) {
            RateShifts r = ((Skygrowth) ne).rateShiftsInput.get();
            // reinitialize the differences array such that all values are 0
            differences = new double[NeFunctionInput.get().neDynamicsInput.get().size()*r.getDimension()];
            // for each time point define in hte rate shifts, get the Ne for that deme
            for (int j = 0; j < r.getDimension(); j++) {
                double time = r.getValue(j);
                // get the corresponding interval in rate shifts
                int interval = 0;
                for (int k = 0; k < intTimes.length; k++) {
                    if (time > intTimes[k]) {
                        interval = k;
                        break;
                    }
                }

                for (int i = 0; i < NeFunctionInput.get().neDynamicsInput.get().size(); i++){
                    // get the covariates for the
                    differences[j*NeFunctionInput.get().neDynamicsInput.get().size()+i] =
                            getRates(j, i) -
                                    Math.log(NeFunctionInput.get().neDynamicsInput.get().get(i).getNeInterval(j));
                }
            }
            // standardize the differences such that the mean is 0 and then compute the logP
            mean = 0;
            for (int i = 0; i < differences.length; i++) {
                mean += differences[i];
            }
            mean/=differences.length;
            for (int i = 0; i < differences.length; i++) {
                differences[i] -= mean;
                logP += dist.logDensity(differences[i]);
            }
        }
        return logP;
    }

    private double getRates(int i, int k) {
        double lograte = 0;
        for (int j = 0; j < covariateListInput.get().size(); j++) {
            if (indicatorInput.get().getArrayValue(j) > 0.0) {
                lograte += scalerInput.get().getArrayValue(j)
                        * covariateListInput.get().get(j).getArrayValue(verticalEntries * i + k);
            }
        }
        return lograte;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {

    }
}
