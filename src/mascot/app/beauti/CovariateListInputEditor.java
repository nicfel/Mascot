package mascot.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.BooleanParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.collections.FXCollections;
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

//	private ObservableList<CovariateRow> covariateObList;

	TableView<CovariateRow> table;

	CovariateList covariateList;

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
		pane.getChildren().clear();
		pane = FXUtils.newHBox();

		ObservableList<CovariateRow> covariateObList = FXCollections.observableArrayList();
		covariateList = (CovariateList) input.get();

//TODO not used. What they do?
//		JDialog dialog = new JDialog((JDialog) null, true);
//		dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
//		dialog.setLocationRelativeTo(null);
//		dialog.setTitle("Isolation with Migration Annotator");

		Button addButton = new Button("Add predictor from file");
		Button updateButton = new Button("Update settings");

		addButton.setOnAction(e -> {
			// TODO cannot set single selection
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
						beastfx.app.util.Alert.showMessageDialog(this,
								error,"Predictor parsing error",
								beastfx.app.util.Alert.ERROR_MESSAGE);
					}
				} catch (IOException ex) {
					ex.printStackTrace();

				}

				// add a transformation vector
				addTransfVector(covariateList);

				try {
					covariateList.initAndValidate();
				} catch (Exception ex) {
					System.err.println("Error initializing covariates.");
					beastfx.app.util.Alert.showMessageDialog(this,
							ex.getMessage(),"Error initializing covariates",
							beastfx.app.util.Alert.ERROR_MESSAGE);
				}

				// data to table
				covariateListToTable(covariateList);
			} // end if file

		});

		VBox boxVert = FXUtils.newVBox();

		HBox boxHoriz = FXUtils.newHBox();
		boxHoriz.getChildren().add(addButton);
		boxHoriz.getChildren().add(updateButton);
//		boxHoriz.add(Box.createHorizontalGlue());
		boxVert.getChildren().add(boxHoriz);

		table = new TableView<>();
		table.setEditable(true);
		table.setPrefWidth(600);
		table.setPrefHeight(150);
		table.getSelectionModel().setSelectionMode(SelectionMode.SINGLE);

		TableColumn<CovariateRow, String> column1 = new TableColumn<>("predictor name");
		column1.setCellValueFactory(new PropertyValueFactory<>("id"));
		column1.setPrefWidth(250);
		column1.setEditable(false);

		TableColumn<CovariateRow, Boolean> column2 = new TableColumn<>("transform");
		column2.setPrefWidth(100);
		column2.setEditable(true);
		column2.setCellValueFactory(new PropertyValueFactory<>("transformed"));
//		column2.setCellFactory(CheckBoxTableCell.forTableColumn(column2));
		column2.setCellFactory(CheckBoxTableCell.forTableColumn(i -> covariateObList.get(i).transformed));

		TableColumn<CovariateRow, String> column3 = new TableColumn<>("time dependent");
		column3.setCellValueFactory(new PropertyValueFactory<>("timeDependent"));
		column3.setPrefWidth(150);
		column3.setEditable(false);//TODO check
//		column3.setCellFactory(TextFieldTableCell.forTableColumn());
//		column3.setOnEditCommit(
//				t -> t.getTableView().getItems().
//						get(t.getTablePosition().getRow()).
//						setTimeDependent(t.getNewValue())
//		);

		TableColumn<CovariateRow, Boolean> column4 = new TableColumn<>("remove predictor");
		column4.setPrefWidth(100);
		column4.setEditable(true);
		column4.setCellValueFactory(new PropertyValueFactory<>("removePredictor"));
//		column4.setCellFactory(CheckBoxTableCell.forTableColumn(column4));
		column4.setCellFactory(CheckBoxTableCell.forTableColumn(i -> covariateObList.get(i).removePredictor));


		table.getColumns().add(column1);
		table.getColumns().add(column2);
		table.getColumns().add(column3);
		table.getColumns().add(column4);

		table.setItems(covariateObList);

		covariateListToTable(covariateList);

		ScrollPane scrollPane = new ScrollPane(table);
//		scrollPane.setPreferredSize(new Dimension(450, 110));

		boxVert.getChildren().add(scrollPane);

		// adds a new dummy variable to the covariates list
		updateButton.setOnAction(e -> {
			// Note: table.getItems() get ObservableList<>, not the GUI selections.
			// This is completely different to Swing table.getModel().getValueAt.
			if (table.getItems().size() != covariateList.covariatesInput.get().size())
				throw new IllegalArgumentException();

			for (int i = covariateList.covariatesInput.get().size() - 1; i >= 0; i--) {
				if (table.getItems().get(i).isRemovePredictor()) {
					covariateList.covariatesInput.get().remove(i);
				} else {
					covariateList.covariatesInput.get().get(i).transformed = table.getItems().get(i).isTransformed();
				}
			}
            // add a transformation vector
			addTransfVector(covariateList);

			try {
				covariateList.initAndValidate();
			} catch (Exception ex) {
				System.err.println("Error initializing covariates.");
				beastfx.app.util.Alert.showMessageDialog(this,
						ex.getMessage(),"Error initializing covariates",
						beastfx.app.util.Alert.ERROR_MESSAGE);
			}

			// data to table
			covariateListToTable(covariateList);
		});


		// add a transformation vector
		addTransfVector(covariateList);

		pane.getChildren().add(boxVert);
		getChildren().add(pane);
	}

	private void addTransfVector(CovariateList covariateList) {
		// add a transformation vector
		Boolean[] transform = new Boolean[covariateList.covariatesInput.get().size()];
		for (int i = 0; i < transform.length; i++)
			transform[i] = covariateList.covariatesInput.get().get(i).transformed;

		covariateList.transformInput.set(new BooleanParameter(transform));
	}

	// push covariateList.covariatesInput into Table
	private void covariateListToTable(CovariateList covariateList) {
		table.getItems().clear();
		// add data
		List<Covariate> covs = covariateList.covariatesInput.get();
		for (Covariate cov : covs) {
			CovariateRow cw = new CovariateRow(cov.getID(), cov.transformed, cov.isTimeDependent);
			table.getItems().add(cw);
		}
		table.refresh();
		refreshPanel();
	}


	/**
	 * table model
	 */
	public class CovariateRow {
		String id;
		BooleanProperty transformed;
		String timeDependent;
		BooleanProperty removePredictor;

		public CovariateRow(String id, boolean transformed, String timeDependent) {
			this.id = id;
			this.transformed = new SimpleBooleanProperty(transformed);
			this.timeDependent = timeDependent;
			this.removePredictor = new SimpleBooleanProperty(false);

			this.transformed.addListener((ov, t, t1) -> {
//					System.out.println(id + " transformed : " + t + " => " + t1);
				setTransformed(t1);
			});

			this.removePredictor.addListener((ov, t, t1) -> {
//					System.out.println(id + " removePredictor : " + t + " => " + t1);
				setRemovePredictor(t1);
			});
		}

		public String getId() {
			return id;
		}

		public void setId(String id) {
			this.id = id;
		}

		public boolean isTransformed() {
			return transformed.get();
		}

		public void setTransformed(boolean transformed) {
			this.transformed.set(transformed);
		}

		public String getTimeDependent() {
			return timeDependent;
		}

		public void setTimeDependent(String timeDependent) {
			this.timeDependent = timeDependent;
		}

		public boolean isRemovePredictor() {
			return removePredictor.get();
		}

		public void setRemovePredictor(boolean removePredictor) {
			this.removePredictor.set(removePredictor);
		}
	}

}
