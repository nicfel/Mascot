package mascot.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.BooleanParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.Alert;
import beastfx.app.util.FXUtils;
import javafx.collections.ObservableList;
import javafx.scene.control.*;
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import mascot.glmmodel.Covariate;
import mascot.glmmodel.CovariateList;
import mascot.util.FolderUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class CovariateListInputEditor extends InputEditor.Base {
//	CovariateList covariateList;
	private ObservableList<CovariateWrapper> covariateObList;

	TableView<CovariateWrapper> table;

	public CovariateListInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return CovariateList.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption bExpandOption,
			boolean addButtons) {
		super.init(input, beastObject, itemNr, ExpandOption.TRUE, addButtons);

		CovariateList covariateList = (CovariateList) input.get();

//TODO not used. What they do?
//		JDialog dialog = new JDialog((JDialog) null, true);
//		dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
//		dialog.setLocationRelativeTo(null);
//		dialog.setTitle("Isolation with Migration Annotator");

		Button addButton = new Button("Add predictor from file");
		Button updateButton = new Button("Update settings");

		addButton.setOnAction(e -> {
			File file = FXUtils.getLoadFile("Select predictor file", FolderUtils.getLastDir(), null, "csv");
			if (file != null) {

				FolderUtils.setLastDirToParentOf(file);
				
				// Read in File
				BufferedReader reader;
				String error = "";
				try {
					reader = new BufferedReader(new FileReader(file));
					String line = reader.readLine();
					List<String> rawValues = new ArrayList<>();
					while (line != null) {
						if (line.contains(","))
							rawValues.add(line);
						// read next line
						line = reader.readLine();
					}
					reader.close();

					// get the id of the predictor
				    String id = file.getName().replace(".csv", "");
					Covariate newCov = new Covariate(rawValues, id);
					covariateList.covariatesInput.get().add(newCov);
					// String matching to see if it is a migration or Ne predictor
					if (covariateList.getID().contains("migration"))
						error = covariateList.initMigrationFromRawValues(covariateList.covariatesInput.get().size()-1);
					else
						error = covariateList.initNeFromRawValues(covariateList.covariatesInput.get().size()-1);
					
					if (!error.contentEquals("")) {
						Alert.showMessageDialog(this,
								error,"Predictor parsing error",
								Alert.ERROR_MESSAGE);
					}
				} catch (IOException ex) {
					ex.printStackTrace();

				}
			}
			try {
				covariateList.initAndValidate();
			} catch (Exception ex) {
				System.err.println("Error initializing covariates.");
			}

			//TODO update table

			refreshPanel();
		});

//		String[] columnNames = { "predictor name", "transform", "time dependent", "remove predictor" };
//		Object[][] data = new Object[covariateList.covariatesInput.get().size()][4];
//		for (int i = 0; i < covariateList.covariatesInput.get().size(); i++) {
//			data[i][0] = covariateList.covariatesInput.get().get(i).getID();
//			data[i][1] = covariateList.covariatesInput.get().get(i).transformed;
//			data[i][2] = covariateList.covariatesInput.get().get(i).isTimeDependent;
//			data[i][3] = false;
//		}

		VBox boxVert = FXUtils.newVBox();

		HBox boxHoriz = FXUtils.newHBox();
		boxHoriz.getChildren().add(addButton);
		boxHoriz.getChildren().add(updateButton);
//		boxHoriz.add(Box.createHorizontalGlue());
		boxVert.getChildren().add(boxHoriz);

//		JTable table = new JTable(data, columnNames);
		table = new TableView<>();
		table.setEditable(true);
		table.setPrefWidth(1024);
		table.getSelectionModel().setSelectionMode(SelectionMode.MULTIPLE);

		TableColumn<CovariateWrapper, String> column1 = new TableColumn<>("predictor name");
		column1.setCellValueFactory(new PropertyValueFactory<>("name"));
		column1.setPrefWidth(250);
		column1.setEditable(false);

		TableColumn<CovariateWrapper, Boolean> column2 = new TableColumn<>("transform");
		column2.setPrefWidth(250);
		column2.setEditable(true);
		column2.setCellValueFactory(new PropertyValueFactory<>("time dependent"));
		column2.setCellFactory(CheckBoxTableCell.forTableColumn(column2));

		TableColumn<CovariateWrapper, String> column3 = new TableColumn<>("time dependent");
		column3.setCellValueFactory(new PropertyValueFactory<>("name"));
		column3.setPrefWidth(250);
		column3.setEditable(true);

		TableColumn<CovariateWrapper, Boolean> column4 = new TableColumn<>("remove predictor");
		column4.setPrefWidth(250);
		column4.setEditable(true);
		column4.setCellValueFactory(new PropertyValueFactory<>("location"));
		column4.setCellFactory(CheckBoxTableCell.forTableColumn(column4));


//		table.getColumn("transform").setCellRenderer(new CheckBoxRenderer());
//		table.getColumn("transform").setCellEditor(new DefaultCellEditor(new JCheckBox()));
//		table.getColumn("remove predictor").setCellRenderer(new CheckBoxRenderer());
//		table.getColumn("remove predictor").setCellEditor(new DefaultCellEditor(new JCheckBox()));

		List<Covariate> covs = covariateList.covariatesInput.get();
		for (Covariate cov : covs) {
			CovariateWrapper cw = new CovariateWrapper(cov.getID(), cov.transformed, cov.isTimeDependent, false);
			covariateObList.add(cw);
		}
		table.setItems(covariateObList);

		ScrollPane scrollPane = new ScrollPane(table);
//		scrollPane.setPreferredSize(new Dimension(450, 110));

		boxVert.getChildren().add(scrollPane);

		// adds a new dummy variable to the covariates list
		updateButton.setOnAction(e -> {
			assert covariateObList.size() == covariateList.covariatesInput.get().size();

			for (int i = covariateList.covariatesInput.get().size() - 1; i >= 0; i--) {
				if (covariateObList.get(i).isRemovePredictor()) {
					covariateList.covariatesInput.get().remove(i);
				} else {
					covariateList.covariatesInput.get().get(i).transformed = covariateObList.get(i).isTransformed();
				}
			}

			// add a transformation vector
			Boolean[] transform = new Boolean[covariateList.covariatesInput.get().size()];
			for (int i = 0; i < transform.length; i++)
				transform[i] = covariateList.covariatesInput.get().get(i).transformed;
			
			covariateList.transformInput.set(new BooleanParameter(transform));

			assert covariateObList.size() == covariateList.covariatesInput.get().size();
			
			//TODO update table
			refreshPanel();
		});
		
		
		// add a transformation vector
		Boolean[] transform = new Boolean[covariateList.covariatesInput.get().size()];
		for (int i = 0; i < transform.length; i++)
			transform[i] = covariateList.covariatesInput.get().get(i).transformed;
		
		covariateList.transformInput.set(new BooleanParameter(transform));

		this.pane = FXUtils.newHBox();
		getChildren().add(pane);

		pane.getChildren().add(boxVert);

	}

	/**
	 * table model
	 */
	public class CovariateWrapper {
		String id;
		boolean transformed;
		String timeDependent;
		boolean removePredictor;

		public CovariateWrapper(String id, boolean transformed, String timeDependent, boolean removePredictor) {
			this.id = id;
			this.transformed = transformed;
			this.timeDependent = timeDependent;
			this.removePredictor = removePredictor;
		}

		public String getId() {
			return id;
		}

		public void setId(String id) {
			this.id = id;
		}

		public boolean isTransformed() {
			return transformed;
		}

		public void setTransformed(boolean transformed) {
			this.transformed = transformed;
		}

		public String getTimeDependent() {
			return timeDependent;
		}

		public void setTimeDependent(String timeDependent) {
			this.timeDependent = timeDependent;
		}

		public boolean isRemovePredictor() {
			return removePredictor;
		}

		public void setRemovePredictor(boolean removePredictor) {
			this.removePredictor = removePredictor;
		}
	}

/*
	public class CheckBoxRenderer extends JCheckBox implements TableCellRenderer {

		CheckBoxRenderer() {
		}

		public Component getTableCellRendererComponent(JTable table, Object value,
            boolean isSelected, boolean hasFocus, int row, int column) {        	
			setSelected((value != null && ((Boolean) value).booleanValue()));
			return this;
        }
	}

	class CheckEditor extends DefaultCellEditor {
		protected JCheckBox button;

		private String label;

		private boolean value;

		public CheckEditor(JCheckBox checkBox, boolean value) {
			super(checkBox);
			button = new JCheckBox();
			button.setText("Enable logging");
			button.setSelected(false);
			this.value = value;

			button.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					fireEditingStopped();
				}
			});
		}

		public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row,
				int column) {
			if (isSelected) {
				setForeground(table.getSelectionForeground());
				// super.setBackground(table.getSelectionBackground());
				setBackground(table.getSelectionBackground());
			} else {
				setForeground(table.getForeground());
				setBackground(table.getBackground());
			}
			this.value = !this.value;
			return button;
		}

		public Object getCellEditorValue() {
			return value;
		}

		public boolean stopCellEditing() {
			return super.stopCellEditing();
		}

		protected void fireEditingStopped() {
			super.fireEditingStopped();
		}
	}
*/
}
