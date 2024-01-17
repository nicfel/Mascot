package mascot.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.operator.AdaptableOperatorSampler;
import beast.base.evolution.operator.kernel.BactrianScaleOperator;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.distribution.Prior;
import beast.base.inference.distribution.Uniform;
import beast.base.inference.operator.kernel.BactrianRandomWalkOperator;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiSubTemplate;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.control.*;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import mascot.dynamics.RateShifts;
import mascot.parameterdynamics.*;
import mascot.util.Difference;
import mascot.util.First;

import java.util.ArrayList;
import java.util.List;

public class NeDynamicsListInputEditor extends InputEditor.Base {

	NeDynamicsList neDynamicsList;
	
    public BeautiSubTemplate hyperPriorTemplate = null;

	
	public NeDynamicsListInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return NeDynamicsList.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption,
			boolean bAddButtons) {
		
		
		
		neDynamicsList = (NeDynamicsList) input.get();
		
		String pId = neDynamicsList.getID().substring(neDynamicsList.getID().indexOf(".t:")+3, neDynamicsList.getID().length());

		
		VBox boxVert = FXUtils.newVBox();
		for (int i = 0; i < neDynamicsList.neDynamicsInput.get().size(); i++) {
			
			Label dynlablel = new Label(neDynamicsList.neDynamicsInput.get().get(i).getID().replace("NeDynamics.", "Ne Dynamics in ").replace(".t:" + pId, ""));
			String state = neDynamicsList.neDynamicsInput.get().get(i).getID();
			
			ObservableList<String> skylineExamples = FXCollections.observableArrayList();
			
			skylineExamples.addAll(List.of(
					"NotSet",
					"Skyline",
					"Constant",
					"Exponential"));

			
	        ComboBox<String> dynamicsComboBox = new ComboBox<>(skylineExamples);
	        dynamicsComboBox.setTooltip(new Tooltip("Set the Ne Dynamics for this state"));
	        dynamicsComboBox.setEditable(true);

	        
        	if(neDynamicsList.neDynamicsInput.get().get(i) instanceof NotSet) {	   
        		dynamicsComboBox.getSelectionModel().select(0);
        	}else if (neDynamicsList.neDynamicsInput.get().get(i) instanceof Skygrowth) {
        		dynamicsComboBox.getSelectionModel().select(1);
        	}else if (neDynamicsList.neDynamicsInput.get().get(i) instanceof ConstantNe) {
        		dynamicsComboBox.getSelectionModel().select(2);
        	}else if (neDynamicsList.neDynamicsInput.get().get(i) instanceof ExponentialNe) {
        		dynamicsComboBox.getSelectionModel().select(3);
        	}

        	
			
        	dynamicsComboBox.setOnAction(e -> {				
	        	switch (dynamicsComboBox.getSelectionModel().getSelectedItem()) {
	    			case "NotSet" : 
	    				setToNotSet(state);
	    				break;
	        		case "Skyline" : 
	        			setToSkyline(state);
	        			break;
	        		case "Constant" :
	        			setToConstant(state);
	        			break;
	        		case "Exponential" :
	        			setToExponential(state);
	        			break;
	        	}
	        	refreshPanel();
	        	sync();
	        });    
			
        	
			
			HBox boxHoriz = FXUtils.newHBox();
			boxHoriz.getChildren().add(dynlablel);
			if (dynlablel.toString().replace("Ne Dynamics in ", "").contains("NOT_SET")) {
				Label label = new Label("GUESS THE TYPES OF THE SAMPLES FIRST");
				boxHoriz.getChildren().add(label);

			}else {
				boxHoriz.getChildren().add(dynamicsComboBox);
			}
			
    		
			
			ObservableList<Integer> nrRateShiftsExamples = FXCollections.observableArrayList();
			
			
			for (int r = 2; r < 100; r++)
				nrRateShiftsExamples.add(r);

			
	        ComboBox<Integer> rateShiftsComboBox = new ComboBox<>(nrRateShiftsExamples);
	        rateShiftsComboBox.setTooltip(new Tooltip("Set the number of rate shifts for the skygrowth dynamics"));
	        rateShiftsComboBox.setEditable(true);
			
        	if (neDynamicsList.neDynamicsInput.get().get(i) instanceof Skygrowth) {
        		int nrShifts = ((Skygrowth) neDynamicsList.neDynamicsInput.get().get(i)).rateShiftsInput.get().valuesInput.get().size();
        		
        		rateShiftsComboBox.getSelectionModel().select(nrShifts-2);        		
        		Label intervalLabel = new Label("number of Ne's to estimate");
        		
        		rateShiftsComboBox.setOnAction(e -> {
    				
        			int nrIntervals = rateShiftsComboBox.getSelectionModel().getSelectedItem();
        			setRateShifts(state, nrIntervals);    
        			    	        	
        			refreshPanel();
    	        });        		
        		boxHoriz.getChildren().add(intervalLabel);
        		boxHoriz.getChildren().add(rateShiftsComboBox);
        	}
			boxVert.getChildren().add(boxHoriz);		
			
			
		}

		this.pane = FXUtils.newHBox();
		getChildren().add(pane);

		pane.getChildren().add(boxVert);	
		sync();
		refreshPanel();
	}

	private void setToNotSet(String state) {
		MCMC mcmc = (MCMC) doc.mcmc.get();

		int currentIndex=-1;
		for (int i = 0; i < neDynamicsList.neDynamicsInput.get().size(); i++) {
			if (neDynamicsList.neDynamicsInput.get().get(i).getID().contentEquals(state))
				currentIndex = i;
		}	
		
		removeParameters(neDynamicsList.neDynamicsInput.get().get(currentIndex), mcmc);

		String id = neDynamicsList.neDynamicsInput.get().get(currentIndex).getID();
		
    	NotSet notSet = new NotSet();
    	notSet.setID(id);
    	neDynamicsList.neDynamicsInput.get().set(currentIndex, notSet);		
	}
	
	private void setToSkyline(String state) {
		MCMC mcmc = (MCMC) doc.mcmc.get();

		int currentIndex=-1;
		for (int i = 0; i < neDynamicsList.neDynamicsInput.get().size(); i++) {
			if (neDynamicsList.neDynamicsInput.get().get(i).getID().contentEquals(state))
				currentIndex = i;
		}

		removeParameters(neDynamicsList.neDynamicsInput.get().get(currentIndex), mcmc);

		String id = neDynamicsList.neDynamicsInput.get().get(currentIndex).getID();		
		String pId = id.substring(id.indexOf(".t:")+3, id.length());
		String location = id.substring(id.indexOf(".")+1, id.indexOf(".t:"));
		
		Skygrowth skygrowth = new Skygrowth();
		skygrowth.setID(id);
		
		RateShifts rateShifts = new RateShifts();
		
		
		String[] tmp = state.split("\\.");
		
		rateShifts.setID("SkygrowthRateShifts." + tmp[1]);
		List<Double> vals = new ArrayList<>();
//		vals.add(0.5);
		vals.add(1.0);
		
		rateShifts.initByName("value", vals, "tree", doc.pluginmap.get("Tree.t:"+pId));
		
		RealParameter logNe = new RealParameter();		
		logNe.setID("SkylineNe." + location + ".t:"+pId);
		logNe.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "0.0", 1);
		
		addParameter(logNe, pId, "SkylineNe."+location, mcmc, true);	
		
		skygrowth.initByName("rateShifts", rateShifts, "logNe", logNe);          				
    	neDynamicsList.neDynamicsInput.get().set(currentIndex, skygrowth);

    	
	}
	
	private void setToConstant(String state) {
		MCMC mcmc = (MCMC) doc.mcmc.get();

		
		int currentIndex=-1;
		for (int i = 0; i < neDynamicsList.neDynamicsInput.get().size(); i++) {
			if (neDynamicsList.neDynamicsInput.get().get(i).getID().contentEquals(state))
				currentIndex = i;
		}		
		
		removeParameters(neDynamicsList.neDynamicsInput.get().get(currentIndex), mcmc);

		String id = neDynamicsList.neDynamicsInput.get().get(currentIndex).getID();		
		String pId = id.substring(id.indexOf(".t:")+3, id.length());
		String location = id.substring(id.indexOf(".")+1, id.indexOf(".t:"));
		
		ConstantNe constant = new ConstantNe();
		constant.setID(id);
		
		RealParameter Ne = new RealParameter();		
		Ne.setID("ConstantNe." + location + ".t:"+pId);
		Ne.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "0.0", 1);
		
		addParameter(Ne, pId, "ConstantNe."+location, mcmc, false);	
		
		constant.initByName("logNe", Ne);          				
    	neDynamicsList.neDynamicsInput.get().set(currentIndex, constant);
	}

	private void setToExponential(String state) {
		MCMC mcmc = (MCMC) doc.mcmc.get();
		
		int currentIndex=-1;
		for (int i = 0; i < neDynamicsList.neDynamicsInput.get().size(); i++) {
			if (neDynamicsList.neDynamicsInput.get().get(i).getID().contentEquals(state))
				currentIndex = i;
		}
		removeParameters(neDynamicsList.neDynamicsInput.get().get(currentIndex), mcmc);

		
		String id = neDynamicsList.neDynamicsInput.get().get(currentIndex).getID();
		String pId = id.substring(id.indexOf(".t:")+3, id.length());
		String location = id.substring(id.indexOf(".")+1, id.indexOf(".t:"));		
	
		ExponentialNe exponential = new ExponentialNe();
		exponential.setID(id);
		
		RealParameter NeNull = new RealParameter();		
		NeNull.setID("NeNull." +  location + ".t:"+pId);
		NeNull.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "0.0", 1);
		
		RealParameter GrowthRate = new RealParameter();		
		GrowthRate.setID("GrowthRate."  + location + ".t:"+pId);
		GrowthRate.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "0.0", 1);
		
		addParameter(NeNull, pId, "NeNull."+location, mcmc, false);
		addParameter(GrowthRate, pId, "GrowthRate."+location, mcmc, false);
		
		exponential.initByName("NeNull", NeNull, "growthRate", GrowthRate);
    	neDynamicsList.neDynamicsInput.get().set(currentIndex, exponential);   
	}
	
	private void setToLogistic(String state) {
		MCMC mcmc = (MCMC) doc.mcmc.get();
		
		int currentIndex=-1;
		for (int i = 0; i < neDynamicsList.neDynamicsInput.get().size(); i++) {
			if (neDynamicsList.neDynamicsInput.get().get(i).getID().contentEquals(state))
				currentIndex = i;
		}
		removeParameters(neDynamicsList.neDynamicsInput.get().get(currentIndex), mcmc);

		
		String id = neDynamicsList.neDynamicsInput.get().get(currentIndex).getID();
		String pId = id.substring(id.indexOf(".t:")+3, id.length());
		String location = id.substring(id.indexOf(".")+1, id.indexOf(".t:"));		
	
		LogisticNe logistic = new LogisticNe();
		logistic.setID(id);
		
		RealParameter NeNull = new RealParameter();		
		NeNull.setID("Capacity." +  location + ".t:"+pId);
		NeNull.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "0.0", 1);
		
		RealParameter GrowthRate = new RealParameter();		
		GrowthRate.setID("GrowthRate."  + location + ".t:"+pId);
		GrowthRate.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "0.0", 1);
		
		
		RealParameter CaryingProportion = new RealParameter();		
		CaryingProportion.setID("CarryingProportion." +  location + ".t:"+pId);
		CaryingProportion.init(1e-6, 1.0-1e-6, "0.5", 1);

		
		addParameter(NeNull, pId, "Capacity."+location, mcmc, false);
		addParameter(GrowthRate, pId, "GrowthRate."+location, mcmc, false);
		addParameter(CaryingProportion, pId, "CarryingProportion."+location, mcmc, false);

		logistic.initByName("capacity", NeNull, "growthRate", GrowthRate, "carryingProportion", CaryingProportion);
		
    	neDynamicsList.neDynamicsInput.get().set(currentIndex, logistic);   

	}

	private void setRateShifts(String state, int intervals) {
		
		int currentIndex=-1;
		for (int i = 0; i < neDynamicsList.neDynamicsInput.get().size(); i++) {
			if (neDynamicsList.neDynamicsInput.get().get(i).getID().contentEquals(state))
				currentIndex = i;
		}
		
		double stepsize = 1.0 / (intervals - 1);
		List<Double> vals = new ArrayList<>();
		double val = 0.0;
		while(val<(1-stepsize/2)) {
			val += stepsize;
			vals.add(val);
		}
		vals.add(1.0);
		
		RateShifts rateShifts =  ((Skygrowth) neDynamicsList.neDynamicsInput.get().get(currentIndex)).rateShiftsInput.get();
		rateShifts.valuesInput = new Input("value", "input of timings of rate shifts relative to the most recent sample", vals);
//		rateShifts.initByName("value", vals, "tree", rateShifts.treeInput.get());

	}
	
	private void addParameter(RealParameter parameter, String pId, String dynamics, MCMC mcmc, boolean isSkykline) {
		doc.addPlugin(parameter);
		
		// add the parameter as a new state
		doc.connect(parameter, "state", "stateNode");	
		doc.connect(parameter, "tracelog", "log");
		// add the operator to the amvn		
		doc.connect(parameter, "AVMNNoTransform.Mascot." + pId, "f");
		
		// add to operator
		AdaptableOperatorSampler ops = new AdaptableOperatorSampler();
		ops.setID(dynamics + ".Scaler.t:" + pId);
		ops.paramInput.setValue(parameter, ops);
		ops.m_pWeight.setValue(1.0, ops);
		
//		BactrianScaleOperator bso = new BactrianScaleOperator();
//		bso.setID(dynamics + ".ScalerX.t:" + pId);
//		bso.parameterInput.setValue(parameter, bso);
//		bso.scaleFactorInput.setValue(0.5, bso);
//		bso.scaleUpperLimit.setValue(10.0, bso);
//		bso.m_pWeight.setValue(1.0, bso);

		BactrianRandomWalkOperator bso = new BactrianRandomWalkOperator();
		bso.setID(dynamics + ".ScalerX.t:" + pId);
		bso.parameterInput.setValue(parameter, bso);
		bso.scaleFactorInput.setValue(0.5, bso);
		bso.m_pWeight.setValue(1.0, bso);
		
		// add the bactrian to the adaptable		
		ops.operatorsInput.get().add(bso);
		// add the amvnm to the adaptable		
		ops.operatorsInput.get().add((Operator) doc.pluginmap.get("AVMNOperator.Mascot."+pId));

		doc.addPlugin(bso);
		doc.addPlugin(ops);

		
		// add the adaptable to mcmc
		doc.connect(ops, "mcmc", "operator");

    	Prior paramPrior = new Prior();
		Uniform uni = new Uniform();		
		uni.setID("Uniform.");

		if (!isSkykline) {
			uni.initByName("lower", Double.NEGATIVE_INFINITY, "upper", Double.POSITIVE_INFINITY);		
			paramPrior.setID(dynamics + ".Prior.t:" + pId);

			paramPrior.distInput.setValue(uni, paramPrior);

			paramPrior.m_x.setValue(parameter, paramPrior);
		}else {
			if (parameter.getID().startsWith("CarryingProportion")) {
				uni.initByName("lower", 0.0, "upper", 1.0);		

			}
			
			uni.initByName("lower", Double.NEGATIVE_INFINITY, "upper", Double.POSITIVE_INFINITY);		
			paramPrior.setID(dynamics + ".Prior.t:" + pId);

			paramPrior.distInput.setValue(uni, paramPrior);

			
			Difference diff = new Difference();
			diff.setID("diff." + parameter.getID());
			diff.initByName("arg", parameter);
			paramPrior.m_x.setValue(diff, paramPrior);			
			
			
	    	Prior paramPrior2 = new Prior();
	    	paramPrior2.setID(dynamics + ".FirstPrior.t:" + pId);	    	
	    	First first = new First();
	    	first.setID("first." + parameter.getID());
	    	first.initByName("arg", parameter);
	    	paramPrior2.m_x.setValue(first, paramPrior2);
	    	
			Uniform uni2 = new Uniform();		
			uni2.setID("Uniform.");
			uni2.initByName("lower", Double.NEGATIVE_INFINITY, "upper", Double.POSITIVE_INFINITY);		
			uni2.initByName("lower", Double.NEGATIVE_INFINITY, "upper", Double.POSITIVE_INFINITY);		
			paramPrior2.distInput.setValue(uni2, paramPrior2);



			doc.registerPlugin(paramPrior2);
			doc.connect(paramPrior2, "prior", "distribution");	

		}
		
		doc.registerPlugin(ops);

		doc.registerPlugin(paramPrior);
		doc.connect(paramPrior, "prior", "distribution");	
	}
	
	private void removeParameters(NeDynamics neDynamics, MCMC mcmc) {
		if (neDynamics instanceof ConstantNe) {
			removeParameter(((CompoundDistribution) mcmc.posteriorInput.get()), ((ConstantNe) neDynamics).NeInput.get());
		}else if (neDynamics instanceof ExponentialNe) {
			removeParameter(((CompoundDistribution) mcmc.posteriorInput.get()), ((ExponentialNe) neDynamics).growthRateInput.get());
			removeParameter(((CompoundDistribution) mcmc.posteriorInput.get()), ((ExponentialNe) neDynamics).NeNullInput.get());			
		}else if (neDynamics instanceof Skygrowth) {
			removeParameter(((CompoundDistribution) mcmc.posteriorInput.get()), ((Skygrowth) neDynamics).NeInput.get());
		}else if (neDynamics instanceof LogisticNe) {
			removeParameter(((CompoundDistribution) mcmc.posteriorInput.get()), ((LogisticNe) neDynamics).capacityInput.get());
			removeParameter(((CompoundDistribution) mcmc.posteriorInput.get()), ((LogisticNe) neDynamics).carryingProportionInput.get());
			removeParameter(((CompoundDistribution) mcmc.posteriorInput.get()), ((LogisticNe) neDynamics).growthRateInput.get());
		}
	}

	private void removeParameter(CompoundDistribution compoundDistribution, RealParameter realParameter) {
		
		CompoundDistribution prior = (CompoundDistribution) doc.pluginmap.get("prior");
		prior.initAndValidate();
		String pId = realParameter.getID().substring(realParameter.getID().indexOf(".t:")+3, realParameter.getID().length());
		String baseName = realParameter.getID().replace(".t:" + pId, "");
		
		doc.disconnect(doc.pluginmap.get(baseName + ".Scaler.t:" + pId), "mcmc", "operator");		
		doc.disconnect(realParameter, "AVMNNoTransform.Mascot."+pId, "f");		
		
		doc.disconnect(realParameter, "tracelog", "log");
		doc.disconnect(realParameter, "state", "stateNode");
		

		// set the parameter to not estimated (will be removed by connector, not cleanest solution, but a working solution)
		for (int i=0; i < prior.pDistributions.get().size(); i++) {
			if (prior.pDistributions.get().get(i).getID().contentEquals(baseName + ".Prior.t:" + pId)) {
				if (((Prior) prior.pDistributions.get().get(i)).m_x.get() instanceof RealParameter) {
					RealParameter rp = ((RealParameter) ((Prior) prior.pDistributions.get().get(i)).m_x.get());
					rp.isEstimatedInput.setValue(false, rp);
				}else{
					Difference fun = ((Difference) ((Prior) prior.pDistributions.get().get(i)).m_x.get());
					RealParameter rp  = (RealParameter) fun.functionInput.get();
					rp.isEstimatedInput.setValue(false, rp);
				}
			}
			if (prior.pDistributions.get().get(i).getID().contentEquals(baseName + ".FirstPrior.t:" + pId)) {
				First fun = ((First) ((Prior) prior.pDistributions.get().get(i)).m_x.get());
				RealParameter rp  = (RealParameter) fun.functionInput.get();
				rp.isEstimatedInput.setValue(false, rp);
			}
		}
		
		doc.connectModel();		
		doc.setUpActivePlugins();		
		sync();
		refreshPanel();

	}
	
//	private List<Double> updateRatesShifts(String from, String to, String mrsi, String by) {
//    	if (rateShifts.dateTimeFormatInput.get().contentEquals("decimal")) {
//    		double fromval = Double.parseDouble(from);
//    		double toval = Double.parseDouble(to);
//    		double mrsival = Double.parseDouble(mrsi);
//    		
//    		double byval = Double.parseDouble(by);
//
//    		
//    		List<Double> newShifts = new ArrayList<>();
//    		double curr_time = fromval-mrsival;
//    		while (toval>=curr_time) {
//       			newShifts.add(curr_time);
//    			curr_time += byval;
//    		}
//    		
//    		newShifts.remove(0);
//    		return newShifts;
//    	}else {
//            DateTimeFormatter formatter = DateTimeFormatter.ofPattern(rateShifts.dateTimeFormatInput.get());
//
//            LocalDate fromval = LocalDate.parse(from, formatter);
//            LocalDate toval = LocalDate.parse(to, formatter);
//    		LocalDate mrsival = LocalDate.parse(mrsi, formatter);
//    		LocalDate curr_date = LocalDate.parse(from, formatter);
//   		
//    		// build Duration
//    		String[] dur_string; 
//    		String[] by_string; 
//   		
//    		int year=0, month=0, day=0;
//    		
//    		if (rateShifts.dateTimeFormatInput.get().contains("-")) {
//    			dur_string = rateShifts.dateTimeFormatInput.get().split("-");
//    			by_string = by.split("-");
//    		}else {
//    			dur_string = rateShifts.dateTimeFormatInput.get().split("/");
//    			by_string = by.split("/");
//    		}
//    		
//    		for (int i = 0; i < dur_string.length; i++) {
//    			if (dur_string[i].contentEquals("yyyy"))
//    				year = Integer.parseInt(by_string[i]);
//				else if (dur_string[i].contentEquals("M"))
//    				month = Integer.parseInt(by_string[i]);
//				else if (dur_string[i].contentEquals("dd"))
//    				day = Integer.parseInt(by_string[i]);
//    		}
//
//    		Period temp = Period.of(year, month, day);
//    			
//
//    		List<Double> newShifts = new ArrayList<>();
//    		double curr_time = -(getYear(fromval)- getYear(mrsival));
//    		double end_time = -(getYear(toval)- getYear(mrsival));
//    		
//    		while (curr_time<=end_time) {
//        		newShifts.add(curr_time);
//    			curr_date = curr_date.minus(temp);	    			
//    			curr_time = -(getYear(curr_date)- getYear(mrsival));
//    		}
//    		// first entry is not a rate shift and is therefore removed
////    		newShifts.remove(0);
//    		return newShifts;
//    	}
//
//	}
//	
//	private Double getYear(LocalDate date) {        
//        return date.getYear() + (date.getDayOfYear()-1.0) / (date.isLeapYear() ? 366.0 : 365.0);
//	}

}