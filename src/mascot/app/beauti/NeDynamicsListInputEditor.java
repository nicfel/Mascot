package mascot.app.beauti;


import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.ListInputEditor;
import mascot.parameterdynamics.EffectivePopulationSizeDynamics;

import java.util.List;



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


    /**TODO
     * add components to box that are specific for the beastObject.
     * By default, this just inserts a label with the beastObject ID
     *
     * @param itemBox box to add components to
     * @param beastObject  beastObject to add

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
    }	*/

    @Override
    protected void addItem() {
        super.addItem();
        sync();
        refreshPanel();
    } // addItem

}
