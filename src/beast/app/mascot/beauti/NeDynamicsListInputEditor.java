package beast.app.mascot.beauti;



import java.awt.Component;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Box;
import javax.swing.JButton;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.app.draw.ListInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.mascot.parameterDynamics.EffectivePopulationSizeDynamics;



public class NeDynamicsListInputEditor extends ListInputEditor {
    private static final long serialVersionUID = 1L;

//    List<JButton> rangeButtons;
//
//    List<JButton> taxonButtons;

	public NeDynamicsListInputEditor(BeautiDoc doc) {
		super(doc);
	}

    @Override
    public Class<?> type() {
        return List.class;
    }

    @Override
    public Class<?> baseType() {
        return EffectivePopulationSizeDynamics.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
    	List<?> list = (List<?>) input.get();        
        //m_buttonStatus = ButtonStatus.NONE;
        super.init(input, beastObject, itemNr, isExpandOption, addButtons);
    }


    /**
     * add components to box that are specific for the beastObject.
     * By default, this just inserts a label with the beastObject ID
     *
     * @param itemBox box to add components to
     * @param beastObject  beastObject to add
     */
    @Override
    protected InputEditor addPluginItem(Box itemBox, BEASTInterface beastObject) {
		try {
	    	int listItemNr = ((List<?>) m_input.get()).indexOf(beastObject);
	    	InputEditor editor = doc.getInputEditorFactory().createInputEditor(m_input, listItemNr, beastObject, false, ExpandOption.FALSE, ButtonStatus.NONE, null, doc);
	    	itemBox.add((Component) editor);
	    	return editor;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return this;
    }	

    @Override
    protected void addItem() {
        super.addItem();
        sync();
        refreshPanel();
    } // addItem

}
