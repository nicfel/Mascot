package mascot.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import javafx.scene.control.Button;
import javafx.scene.layout.HBox;
import mascot.dynamics.GLM;

public class GLMInputEditor extends BEASTObjectInputEditor {

	GLM dyn;
	
	public GLMInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return GLM.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption bExpandOption,
			boolean addButtons) {
		super.init(input, beastObject, itemNr, ExpandOption.TRUE, addButtons);

		dyn = (GLM) input.get();

		Button updateButton = new Button("Update settings");

		// adds a new dummy variable to the covariates list
		updateButton.setOnAction(e -> {
			dyn.initAndValidate();
			super.init(input, beastObject, itemNr, bExpandOption, addButtons);
			refreshPanel();
		});

//		VBox boxVert = FXUtils.newVBox();

		HBox boxHoriz = FXUtils.newHBox();
//        boxHoriz.add(Box.createHorizontalGlue());
		boxHoriz.getChildren().add(updateButton);
//		boxVert.getChildren().add(boxHoriz);

		this.pane = FXUtils.newHBox();
		getChildren().add(pane);

		pane.getChildren().add(boxHoriz);
	}
}
