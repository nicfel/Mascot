package beast.app.mascot.beauti;

import java.time.LocalDate;
import java.time.Period;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JTextField;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.mascot.dynamics.RateShifts;

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
    	rateShifts.dateTimeFormatInput.setValue("decimal", rateShifts);

		JButton updateButton = new JButton("Compute rate shifts");
		
		JLabel FormatLabel = new JLabel("Date format");

        String[] dateFormatExamples = {
        		"decimal",
                "dd/M/yyyy",
                "M/dd/yyyy",
                "yyyy/M/dd",
                "dd-M-yyyy",
                "M-dd-yyyy",
                "yyyy-M-dd"};
        
        JComboBox dateFormatComboBox = new JComboBox<>(dateFormatExamples);
        dateFormatComboBox.setToolTipText("Set format used to parse date values");
        dateFormatComboBox.setEditable(true);
        dateFormatComboBox.setMaximumSize(dateFormatComboBox.getPreferredSize());
        dateFormatComboBox.setSelectedItem(dateFormatComboBox.getSelectedItem());        
        
        dateFormatComboBox.addActionListener(e -> {
        	rateShifts.dateTimeFormatInput.setValue(dateFormatComboBox.getSelectedItem(), rateShifts);
        	System.out.println(rateShifts.dateTimeFormatInput.get());
        });    
		
		JLabel mrsiLabel = new JLabel("Most Recent Sample");
		JLabel fromLabel = new JLabel("From");
		JLabel toLabel = new JLabel("To (further in the past)");
		JLabel intervalLabel = new JLabel("interval length");
		JLabel ratesShiftsLabel = new JLabel("rateShifts");

		JTextField mrsiField = new JTextField(10);
		JTextField fromField = new JTextField(10);
		JTextField toField = new JTextField(10);
		JTextField intervalField = new JTextField(5);
		
		mrsiField.setEditable(true);
		fromField.setEditable(true);
		toField.setEditable(true);
		intervalField.setEditable(true);		
		
		JTextField timeIntervals = new JTextField(30);
		if (rateShifts.valuesInput.get()==null)
			timeIntervals.setText("not computed yet");
		else
			timeIntervals.setText(rateShifts.valuesInput.get().toString());

		mrsiField.setText("");
		fromField.setText("");
		toField.setText("");
		intervalField.setText("");
		
		
		// adds a new dummy variable to the covariates list
		updateButton.addActionListener(e -> {
    		rateShifts.valuesInput.setValue(updateRatesShifts(fromField.getText(), toField.getText(), mrsiField.getText(), intervalField.getText()), rateShifts);
			refreshPanel();
		});

		Box boxVert = Box.createVerticalBox();

		Box boxHoriz = Box.createHorizontalBox();
		boxHoriz.add(Box.createHorizontalGlue());
		boxHoriz.add(mrsiLabel);
		boxHoriz.add(mrsiField);
		boxHoriz.add(fromLabel);
		boxHoriz.add(fromField);
		boxHoriz.add(toLabel);
		boxHoriz.add(toField);
		boxHoriz.add(intervalLabel);
		boxHoriz.add(intervalField);
		
		Box boxHoriz2 = Box.createHorizontalBox();
		boxHoriz2.add(Box.createHorizontalGlue());
		boxHoriz2.add(ratesShiftsLabel);
		boxHoriz2.add(timeIntervals);
		
		Box boxHoriz3 = Box.createHorizontalBox();
		boxHoriz3.add(Box.createHorizontalGlue());
		boxHoriz3.add(updateButton);
		boxHoriz3.add(dateFormatComboBox);

		boxVert.add(boxHoriz3);
		boxVert.add(boxHoriz);
		boxVert.add(boxHoriz2);		
		
		add(boxVert);

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
	
	private Double getYear(LocalDate date) {        
        return date.getYear() + (date.getDayOfYear()-1.0) / (date.isLeapYear() ? 366.0 : 365.0);
	}

}
