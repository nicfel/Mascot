package mascot.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.parser.PartitionContext;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiSubTemplate;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.layout.HBox;
import mascot.parameterdynamics.EffectivePopulationSizeDynamics;

import java.util.List;

public class NeDynamicsInputEditor extends InputEditor.Base {

	public NeDynamicsInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return EffectivePopulationSizeDynamics.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int listItemNr, ExpandOption isExpandOption, boolean addButtons) {
        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr= listItemNr;

        HBox itemBox = FXUtils.newHBox();

        EffectivePopulationSizeDynamics neDyn = (EffectivePopulationSizeDynamics) beastObject;
        String text = neDyn.getID();
        Label label = new Label(text);
//        Font font = label.getFont();
//        Dimension size = new Dimension(font.getSize() * 200 / 13, font.getSize() * 25/13);
//        label.setMinimumSize(size);
//        label.setPreferredSize(size);
        itemBox.getChildren().add(label);

        ObservableList<BeautiSubTemplate> obsBOList = FXCollections.observableArrayList();
        obsBOList.addAll(doc.getInputEditorFactory().getAvailableTemplates(input, neDyn, null, doc));
        ComboBox<BeautiSubTemplate> comboBox = new ComboBox<BeautiSubTemplate>(obsBOList);
        comboBox.setId(text+".distr");

        String id = neDyn.getID();
        //Log.warning.println("id=" + id);
        id = id.substring(0, id.indexOf('.'));
        for (BeautiSubTemplate template : obsBOList) {
            if (template.classInput.get() != null && template.getShortClassName().equals(id)) {
                comboBox.setValue(template); // setSelectedItem
            }
        }
        
        comboBox.setOnAction(e -> {
            @SuppressWarnings("unchecked")
			ComboBox<BeautiSubTemplate> comboBox1 = (ComboBox<BeautiSubTemplate>) e.getSource();

            List<BEASTInterface> list = (List<BEASTInterface>) m_input.get();

            BeautiSubTemplate template = (BeautiSubTemplate) comboBox1.getSelectionModel().getSelectedItem();
            //String id = ((BEASTObject) list.get(item)).getID();
            //String partition = BeautiDoc.parsePartition(id);
//            PartitionContext context = doc.getContextFor((BEASTInterface) list.get(itemNr));
            String newContextName = neDyn.getID().substring(neDyn.getID().indexOf(".t:") + 3, neDyn.getID().length());

//            String newContextName = neDyn.getID().substring(neDyn.getID().indexOf(".", neDyn.getID().indexOf(".t:") + 3)+1, neDyn.getID().length());
            System.out.println(newContextName);
            PartitionContext context = new PartitionContext(newContextName);
            System.out.println("================================    ===== ==== ==== " + context);
            EffectivePopulationSizeDynamics dyn1 = (EffectivePopulationSizeDynamics) list.get(itemNr);
            try {
                template.createSubNet(context, true);
            } catch (Exception e1) {
                e1.printStackTrace();
            }
            sync();
            refreshPanel();
        });
//        JPanel panel = new JPanel();
//        panel.add(comboBox);
//        panel.setMaximumSize(size);
        itemBox.getChildren().add(comboBox);

        this.pane = FXUtils.newVBox();
        getChildren().add(pane);

        pane.getChildren().add(itemBox);
	}

}
