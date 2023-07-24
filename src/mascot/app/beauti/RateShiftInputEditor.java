package mascot.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.control.*;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import mascot.dynamics.GLM;
import mascot.dynamics.RateShifts;

import java.time.LocalDate;
import java.time.Period;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.ArrayList;
import java.util.List;

public class RateShiftInputEditor extends InputEditor.Base {

	RateShifts rateShifts;
	
	public RateShiftInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return RateShifts.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption,
			boolean bAddButtons) {
		
		rateShifts = (RateShifts) input.get();
		// check if these rate shifts are for a glm model
		boolean isGLM = false;
		
		
		
		if (plugin instanceof GLM) {
	    	rateShifts.dateTimeFormatInput.setValue("decimal", rateShifts);
	
			Button computeButton = new Button("Compute rate shifts");
			
	//		Label FormatLabel = new Label("Date format");//TODO not used?
	
			ObservableList<String> dateFormatExamples = FXCollections.observableArrayList();
			dateFormatExamples.addAll(List.of(
					"decimal",
					"dd/M/yyyy",
					"M/dd/yyyy",
					"yyyy/M/dd",
					"dd-M-yyyy",
					"M-dd-yyyy",
					"yyyy-M-dd"));
	
	        ComboBox<String> dateFormatComboBox = new ComboBox<>(dateFormatExamples);
	        dateFormatComboBox.setTooltip(new Tooltip("Set format used to parse date values"));
	        dateFormatComboBox.setEditable(true);
			dateFormatComboBox.getSelectionModel().selectFirst();
	//        dateFormatComboBox.setMaximumSize(dateFormatComboBox.getPreferredSize());
	        
	        dateFormatComboBox.setOnAction(e -> {
	        	rateShifts.dateTimeFormatInput.setValue(dateFormatComboBox.getSelectionModel().getSelectedItem(), rateShifts);
	        	System.out.println(rateShifts.dateTimeFormatInput.get());
	        });    
			
			Label mrsiLabel = new Label("Most Recent Sample");
			Label fromLabel = new Label("From");
			Label toLabel = new Label("To (further in the past)");
			Label intervalLabel = new Label("interval length");
			Label ratesShiftsLabel = new Label("rateShifts");
	
			TextField mrsiField = new TextField();
			mrsiField.setPrefColumnCount(10);
			TextField fromField = new TextField();
			fromField.setPrefColumnCount(10);
			TextField toField = new TextField();
			toField.setPrefColumnCount(10);
			TextField intervalField = new TextField();
			intervalField.setPrefColumnCount(5);
			
			mrsiField.setEditable(true);
			fromField.setEditable(true);
			toField.setEditable(true);
			intervalField.setEditable(true);		
			
			TextField timeIntervals = new TextField();
			timeIntervals.setPrefColumnCount(30);
			if (rateShifts.valuesInput.get()==null)
				timeIntervals.setText("not computed yet");
			else
				timeIntervals.setText(rateShifts.valuesInput.get().toString());
	
			mrsiField.setText("");
			fromField.setText("");
			toField.setText("");
			intervalField.setText("");
	
			// adds a new dummy variable to the covariates list
			computeButton.setOnAction(e -> {
	    		try {
					// suggest trim spaces
	    			rateShifts.valuesInput = new Input("value", "input of timings of rate shifts relative to the most recent sample", updateRatesShifts(fromField.getText().trim(),
							toField.getText().trim(), mrsiField.getText().trim(), intervalField.getText().trim()));
	//				rateShifts.valuesInput.setValue(updateRatesShifts(fromField.getText().trim(),
	//						toField.getText().trim(), mrsiField.getText().trim(), intervalField.getText().trim()), rateShifts);
				} catch (DateTimeParseException dtpe) {
					beastfx.app.util.Alert.showMessageDialog(this,
							dtpe.getMessage(),"Date parsing error",
							beastfx.app.util.Alert.ERROR_MESSAGE);
				}
				refreshPanel();
			});
	
			VBox boxVert = FXUtils.newVBox();
	
			HBox boxHoriz = FXUtils.newHBox();
	//		boxHoriz.add(Box.createHorizontalGlue());
			boxHoriz.getChildren().add(mrsiLabel);
			boxHoriz.getChildren().add(mrsiField);
			boxHoriz.getChildren().add(fromLabel);
			boxHoriz.getChildren().add(fromField);
			boxHoriz.getChildren().add(toLabel);
			boxHoriz.getChildren().add(toField);
			boxHoriz.getChildren().add(intervalLabel);
			boxHoriz.getChildren().add(intervalField);
	
			HBox boxHoriz2 = FXUtils.newHBox();
	//		boxHoriz2.add(Box.createHorizontalGlue());
			boxHoriz2.getChildren().add(ratesShiftsLabel);
			boxHoriz2.getChildren().add(timeIntervals);
	
			HBox boxHoriz3 = FXUtils.newHBox();
	//		boxHoriz3.add(Box.createHorizontalGlue());
			boxHoriz3.getChildren().add(computeButton);
			boxHoriz3.getChildren().add(dateFormatComboBox);
	
			boxVert.getChildren().add(boxHoriz3);
			boxVert.getChildren().add(boxHoriz);
			boxVert.getChildren().add(boxHoriz2);
	
			this.pane = FXUtils.newHBox();
			getChildren().add(pane);
	
			pane.getChildren().add(boxVert);
		}else {
			VBox boxVert = FXUtils.newVBox();
			HBox boxHoriz = FXUtils.newHBox();

			ObservableList<Integer> nrRateShiftsExamples = FXCollections.observableArrayList();
			
			if (rateShifts.valuesInput.get().size()<2) {
    			String pId = plugin.getID().substring(plugin.getID().indexOf(".t:")+3, plugin.getID().length());  
    			List<Double> vals = new ArrayList<>();
    			vals.add(0.0);
    			vals.add(1.0);

				rateShifts.initByName("value", vals, "tree", doc.pluginmap.get("Tree.t:"+pId));
    			setRateShifts(rateShifts, 100);  
			}
			
			for (int i = 50; i < 1001; i++)
				nrRateShiftsExamples.add(i);
			
	        ComboBox<Integer> rateShiftsComboBox = new ComboBox<>(nrRateShiftsExamples);
	        rateShiftsComboBox.setTooltip(new Tooltip("Set the number of rate shifts to compute the changing parameters in a stepwise constant way"));
	        rateShiftsComboBox.setEditable(true);
	        
	        int nrShifts = rateShifts.valuesInput.get().size()-1;
	        
    		rateShiftsComboBox.getSelectionModel().select(nrShifts-50);
    		
    		Label intervalLabel = new Label("Set the number of change points to recompute Ne's and migration rates");
    		
    		rateShiftsComboBox.setOnAction(e -> {
				
    			int nrIntervals = rateShiftsComboBox.getSelectionModel().getSelectedItem();
    			setRateShifts(rateShifts, nrIntervals);    
    			    	        	
    			refreshPanel();
	        });        		
    		boxHoriz.getChildren().add(intervalLabel);
    		boxHoriz.getChildren().add(rateShiftsComboBox);
        	
			boxVert.getChildren().add(boxHoriz);		
			
			this.pane = FXUtils.newHBox();
			getChildren().add(pane);
	
			pane.getChildren().add(boxVert);
		}
	}
	

	private List<Double> updateRatesShifts(String from, String to, String mrsi, String by) {
    	if (rateShifts.dateTimeFormatInput.get().contentEquals("decimal")) {
    		double fromval = Double.parseDouble(from);
    		double toval = Double.parseDouble(to);
    		double mrsival = Double.parseDouble(mrsi);
    		
    		double byval = Double.parseDouble(by);

    		
    		List<Double> newShifts = new ArrayList<>();
    		double curr_time = fromval-mrsival;
    		while (toval>=curr_time) {
       			newShifts.add(curr_time);
    			curr_time += byval;
    		}
    		
    		newShifts.remove(0);
    		return newShifts;
    	}else {
            DateTimeFormatter formatter = DateTimeFormatter.ofPattern(rateShifts.dateTimeFormatInput.get());

            LocalDate fromval = LocalDate.parse(from, formatter);
            LocalDate toval = LocalDate.parse(to, formatter);
    		LocalDate mrsival = LocalDate.parse(mrsi, formatter);
    		LocalDate curr_date = LocalDate.parse(from, formatter);
   		
    		// build Duration
    		String[] dur_string; 
    		String[] by_string; 
   		
    		int year=0, month=0, day=0;
    		
    		if (rateShifts.dateTimeFormatInput.get().contains("-")) {
    			dur_string = rateShifts.dateTimeFormatInput.get().split("-");
    			by_string = by.split("-");
    		}else {
    			dur_string = rateShifts.dateTimeFormatInput.get().split("/");
    			by_string = by.split("/");
    		}
    		
    		for (int i = 0; i < dur_string.length; i++) {
    			if (dur_string[i].contentEquals("yyyy"))
    				year = Integer.parseInt(by_string[i]);
				else if (dur_string[i].contentEquals("M"))
    				month = Integer.parseInt(by_string[i]);
				else if (dur_string[i].contentEquals("dd"))
    				day = Integer.parseInt(by_string[i]);
    		}

    		Period temp = Period.of(year, month, day);
    			

    		List<Double> newShifts = new ArrayList<>();
    		double curr_time = -(getYear(fromval)- getYear(mrsival));
    		double end_time = -(getYear(toval)- getYear(mrsival));
    		
    		while (curr_time<=end_time) {
        		newShifts.add(curr_time);
    			curr_date = curr_date.minus(temp);	    			
    			curr_time = -(getYear(curr_date)- getYear(mrsival));
    		}
    		// first entry is not a rate shift and is therefore removed
//    		newShifts.remove(0);
    		return newShifts;
    	}
	}
	
	private void setRateShifts(RateShifts rateShifts, int nrIntervals) {		
		
		double stepsize = 1.0 / (nrIntervals);
		List<Double> vals = new ArrayList<>();
		double val = 0.0;
		while(val<(1-stepsize/2)) {
			vals.add(val);
			val += stepsize;
		}
		vals.add(1.1);
		
		rateShifts.valuesInput = new Input("value", "input of timings of rate shifts relative to the most recent sample", vals);
    	
	}

	
	private Double getYear(LocalDate date) {        
        return date.getYear() + (date.getDayOfYear()-1.0) / (date.isLeapYear() ? 366.0 : 365.0);
	}

}