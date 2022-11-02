package mascot.app.beauti;



import javax.swing.Box;
import javax.swing.JButton;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
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
		
		dyn = (GLM) input.get();

		JButton updateButton = new JButton("Update settings");

		
		// adds a new dummy variable to the covariates list
		updateButton.addActionListener(e -> {
			dyn.initAndValidate();
			super.init(input, beastObject, itemNr, bExpandOption, addButtons);
			refreshPanel();
		});

		Box boxVert = Box.createVerticalBox();

		Box boxHoriz = Box.createHorizontalBox();
		boxHoriz.add(Box.createHorizontalGlue());
		boxHoriz.add(updateButton);
		boxVert.add(boxHoriz);
		
		add(boxHoriz);


		super.init(input, beastObject, itemNr, bExpandOption, addButtons);
	}
}
