package mascot.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.Node;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.layout.VBox;
import mascot.dynamics.StructuredSkyline;

public class SkylineInputEditor extends BEASTObjectInputEditor {

	StructuredSkyline dyn;
	
	public SkylineInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return StructuredSkyline.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption bExpandOption,
			boolean addButtons) {
		super.init(input, beastObject, itemNr, ExpandOption.TRUE, addButtons);

		dyn = (StructuredSkyline) input.get();

		Button updateButton = new Button("Update settings");

		// adds a new dummy variable to the covariates list
		updateButton.setOnAction(e -> {
			dyn.initAndValidate();
			super.init(input, beastObject, itemNr, ExpandOption.TRUE, addButtons);
			paintPane(updateButton);
//			refreshPanel();
		});

		paintPane(updateButton);
	}

	private void paintPane(final Button updateButton) {
		// redraw everything
		ObservableList<Node> nodes = FXCollections.observableArrayList(pane.getChildren());
		pane.getChildren().clear();
		pane = FXUtils.newHBox();
		VBox boxVert = FXUtils.newVBox();

		//2 nodes: Label, HBox from super.init
		for (Node node : nodes) {
			if (node instanceof Label label) {
				boxVert.getChildren().add(label);
				boxVert.getChildren().add(updateButton);
				pane.getChildren().add(0, boxVert);
//			} else if (node instanceof VBox vBox) { // hard code
//				paintPane(vBox);
//				pane.getChildren().add(vBox);
			} else
				pane.getChildren().add(node);
		}
		getChildren().add(pane);
	}

//	private void paintPane(VBox vBox) {
		// redraw everything
//		ObservableList<Node> nodes = FXCollections.observableArrayList(vBox.getChildren());
//		pane.getChildren().clear();
//		pane = FXUtils.newHBox();
//		VBox boxVert = FXUtils.newVBox();
//
//		//2 nodes: Label, HBox from super.init
//		for (Node node : nodes) {
//			if (updateButton != null && node instanceof Label label) {
//				boxVert.getChildren().add(label);
//				boxVert.getChildren().add(updateButton);
//				pane.getChildren().add(0, boxVert);
//			} else if (updateButton != null && node instanceof VBox vBox) { // hard code
//				paintPane(vBox, null);
//			} else
//				pane.getChildren().add(node);
//		}
//		getChildren().add(pane);
//	}
}
