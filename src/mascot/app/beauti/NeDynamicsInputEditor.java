package mascot.app.beauti;

import java.awt.Dimension;
import java.awt.Font;
import java.util.Arrays;
import java.util.List;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.text.StyledEditorKit.FontSizeAction;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiSubTemplate;
import beast.base.parser.PartitionContext;
import beastfx.app.inputeditor.BEASTObjectDialog;
import beastfx.app.inputeditor.InputEditor;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.EffectivePopulationSizeDynamics;

public class NeDynamicsInputEditor extends InputEditor.Base {
	private static final long serialVersionUID = 1L;

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
		
        Box itemBox = Box.createHorizontalBox();

        EffectivePopulationSizeDynamics neDyn = (EffectivePopulationSizeDynamics) beastObject;
        String text = neDyn.getID();
        JLabel label = new JLabel(text);
        Font font = label.getFont();
        Dimension size = new Dimension(font.getSize() * 200 / 13, font.getSize() * 25/13);
        label.setMinimumSize(size);
        label.setPreferredSize(size);
        itemBox.add(label);

        List<BeautiSubTemplate> availableBEASTObjects = doc.getInputEditorFactory().getAvailableTemplates(input, neDyn, null, doc);
        JComboBox<BeautiSubTemplate> comboBox = new JComboBox<BeautiSubTemplate>(availableBEASTObjects.toArray(new BeautiSubTemplate[]{}));
        comboBox.setName(text+".distr");

        String id = neDyn.getID();
        //Log.warning.println("id=" + id);
        id = id.substring(0, id.indexOf('.'));
        for (BeautiSubTemplate template : availableBEASTObjects) {
            if (template.classInput.get() != null && template.getShortClassName().equals(id)) {
                comboBox.setSelectedItem(template);
            }
        }
        
        comboBox.addActionListener(e -> {
            @SuppressWarnings("unchecked")
			JComboBox<BeautiSubTemplate> comboBox1 = (JComboBox<BeautiSubTemplate>) e.getSource();

            List<BEASTInterface> list = (List<BEASTInterface>) m_input.get();

            BeautiSubTemplate template = (BeautiSubTemplate) comboBox1.getSelectedItem();
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
        JPanel panel = new JPanel();
        panel.add(comboBox);
        panel.setMaximumSize(size);
        itemBox.add(panel);
        

        add(itemBox);
	}

}
