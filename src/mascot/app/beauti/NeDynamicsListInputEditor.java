package mascot.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
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
import mascot.parameterdynamics.NeDynamics;
import mascot.parameterdynamics.NeDynamicsList;

import java.util.List;

public class NeDynamicsListInputEditor extends InputEditor.Base {

//	private ObservableList<CovariateRow> covariateObList;

	TableView<NeRow> table;

	NeDynamicsList NeList;

	public NeDynamicsListInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return NeDynamicsList.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption bExpandOption,
			boolean addButtons) {
		super.init(input, beastObject, itemNr, ExpandOption.TRUE, addButtons);
		
		
		pane.getChildren().clear();
//		pane = FXUtils.newHBox();

		ObservableList<NeRow> covariateObList = FXCollections.observableArrayList();
		NeList = (NeDynamicsList) input.get();


 //		Button updateButton = new Button("Update settings");

//		addButton.setOnAction(e -> {
//			// TODO cannot set single selection
//			File file = FXUtils.getLoadFile("Select predictor file", FolderUtils.getLastDir(), null, "csv");
//			if (file != null) {
//
//				FolderUtils.setLastDirToParentOf(file);
//
//				// Read in File
//				BufferedReader reader;
//				String error = "";
//				try {
//					reader = new BufferedReader(new FileReader(file));
//					String line = reader.readLine();
//					List<String> rawValues = new ArrayList<>();
//					while (line != null) {
//						if (line.contains(","))
//							rawValues.add(line);
//						// read next line
//						line = reader.readLine();
//					}
//					reader.close();
//
//					// get the id of the predictor
//				    String id = file.getName().replace(".csv", "");
//					NeDynamics newCov = new Covariate(rawValues, id);
//					covariateList.covariatesInput.get().add(newCov);
//					// String matching to see if it is a migration or Ne predictor
//					if (covariateList.getID().contains("migration"))
//						error = covariateList.initMigrationFromRawValues(covariateList.covariatesInput.get().size()-1);
//					else
//						error = covariateList.initNeFromRawValues(covariateList.covariatesInput.get().size()-1);
//
//					if (!error.contentEquals("")) {
//						beastfx.app.util.Alert.showMessageDialog(this,
//								error,"Predictor parsing error",
//								beastfx.app.util.Alert.ERROR_MESSAGE);
//					}
//				} catch (IOException ex) {
//					ex.printStackTrace();
//
//				}
//
//				// add a transformation vector
//				addTransfVector(covariateList);
//
//				try {
//					covariateList.initAndValidate();
//				} catch (Exception ex) {
//					System.err.println("Error initializing covariates.");
//					beastfx.app.util.Alert.showMessageDialog(this,
//							ex.getMessage(),"Error initializing covariates",
//							beastfx.app.util.Alert.ERROR_MESSAGE);
//				}
//
//				// data to table
//				covariateListToTable(covariateList);
//			} // end if file
//
//		});

		VBox boxVert = FXUtils.newVBox();

		HBox boxHoriz = FXUtils.newHBox();
//		boxHoriz.getChildren().add(addButton);
//		boxHoriz.getChildren().add(updateButton);
//		boxHoriz.add(Box.createHorizontalGlue());
		boxVert.getChildren().add(boxHoriz);

		table = new TableView<>();
		table.setEditable(true);
		table.setPrefWidth(600);
		table.setPrefHeight(150);
		table.getSelectionModel().setSelectionMode(SelectionMode.SINGLE);

		TableColumn<NeRow, String> column1 = new TableColumn<>("state/location");
		column1.setCellValueFactory(new PropertyValueFactory<>("id"));
		column1.setPrefWidth(250);
		column1.setEditable(false);

		TableColumn<NeRow, Boolean> column2 = new TableColumn<>("dynamics type");
		column2.setPrefWidth(100);
		column2.setEditable(true);
		column2.setCellValueFactory(new PropertyValueFactory<>("type"));
//		column2.setCellFactory(CheckBoxTableCell.forTableColumn(column2));
		column2.setCellFactory(CheckBoxTableCell.forTableColumn(i -> covariateObList.get(i).transformed));


		table.getColumns().add(column1);
		table.getColumns().add(column2);
		table.setItems(covariateObList);

		covariateListToTable(NeList);

		ScrollPane scrollPane = new ScrollPane(table);
//		scrollPane.setPreferredSize(new Dimension(450, 110));

		boxVert.getChildren().add(scrollPane);

		// adds a new dummy variable to the covariates list
//		updateButton.setOnAction(e -> {
//			// Note: table.getItems() get ObservableList<>, not the GUI selections.
//			// This is completely different to Swing table.getModel().getValueAt.
//			if (table.getItems().size() != NeList.neDynamicsInput.get().size())
//				throw new IllegalArgumentException();
//
//			for (int i = NeList.neDynamicsInput.get().size() - 1; i >= 0; i--) {
//				if (table.getItems().get(i).isRemovePredictor()) {
//					NeList.neDynamicsInput.get().remove(i);
//				} else {
//					NeList.neDynamicsInput.get().get(i).transformed = table.getItems().get(i).isTransformed();
//				}
//			}
//
//			try {
//				NeList.initAndValidate();
//			} catch (Exception ex) {
//				System.err.println("Error initializing covariates.");
//				beastfx.app.util.Alert.showMessageDialog(this,
//						ex.getMessage(),"Error initializing covariates",
//						beastfx.app.util.Alert.ERROR_MESSAGE);
//			}
//
//			// data to table
//			covariateListToTable(NeList);
//		});


		pane.getChildren().add(boxVert);
//		getChildren().add(pane);
	}


	// push covariateList.covariatesInput into Table
	private void covariateListToTable(NeDynamicsList neDynamicsList) {
		table.getItems().clear();
		// add data
		List<NeDynamics> nes = neDynamicsList.neDynamicsInput.get();
		for (NeDynamics ne : nes) {
			NeRow cw = new NeRow(ne.getID());
			table.getItems().add(cw);
		}
		table.refresh();
		refreshPanel();
	}


	/**
	 * table model
	 */
	public class NeRow {
		String id;
		BooleanProperty transformed;
		String timeDependent;
		BooleanProperty removePredictor;

		public NeRow(String id) {
			this.id = id;
			this.transformed = new SimpleBooleanProperty(false);
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
