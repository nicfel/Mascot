package beast.app.mascot.beauti;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import javax.swing.Box;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.WindowConstants;
import javax.swing.table.TableCellRenderer;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.mascot.glmmodel.CovariateList;
import beast.mascot.glmmodel.Covariate;

public class CovariateListInputEditor extends InputEditor.Base {
	CovariateList covariateList;
	
	

	List<Covariate> covs;

	public CovariateListInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return CovariateList.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption,
			boolean bAddButtons) {
		covariateList = (CovariateList) input.get();
		covs = covariateList.covariatesInput.get();

		JDialog dialog = new JDialog((JDialog) null, true);
		dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		dialog.setLocationRelativeTo(null);
		dialog.setTitle("Isolation with Migration Annotator");

		JButton addButton = new JButton("Add covariate from file");
		JButton updateButton = new JButton("Update settings");

		JFileChooser inFileChooser = new JFileChooser();
		addButton.addActionListener(e -> {
			inFileChooser.setDialogTitle("Select covariate file");
			inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
			int returnVal = inFileChooser.showOpenDialog(dialog);
			//
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File fname = inFileChooser.getSelectedFile();
				// Read in File
				BufferedReader reader;
				try {
					reader = new BufferedReader(new FileReader(fname));
					String line = reader.readLine();
					while (line != null) {
						String[] names = line.split(":");
						System.out.println(Arrays.toString(names));

						String[] vals_string = names[1].trim().split("\\s+");
						System.out.println(Arrays.toString(vals_string));

						Double[] new_vals = new Double[vals_string.length];
						for (int i = 0; i < new_vals.length; i++) {
							System.out.println(vals_string[i]);
							new_vals[i] = Double.parseDouble(vals_string[i]);
						}

						Covariate newCov = new Covariate();
						covariateList.covariatesInput.get().add(newCov);

						// read next line
						line = reader.readLine();
					}
					reader.close();
				} catch (IOException ex) {
					ex.printStackTrace();
				}
			}
			try {
				covariateList.initAndValidate();
			} catch (Exception ex) {
				System.err.println("Error initializing covariates.");
			}

			refreshPanel();
		});

		String[] columnNames = { "predictor name", "transform", "time invariant", "remove predictor" };

		Object[][] data = new Object[covariateList.covariatesInput.get().size()][4];
		for (int i = 0; i < covariateList.covariatesInput.get().size(); i++) {
			data[i][0] = covariateList.covariatesInput.get().get(i).getID();
			data[i][1] = covariateList.covariatesInput.get().get(i).logTransformInput.get();
			data[i][2] = covariateList.covariatesInput.get().get(i).timeInvariantInput.get();
			data[i][3] = false;
		}

		Box boxVert = Box.createVerticalBox();

		Box boxHoriz = Box.createHorizontalBox();
		boxHoriz.add(Box.createHorizontalGlue());
		boxHoriz.add(addButton);
		boxHoriz.add(updateButton);
		boxVert.add(boxHoriz);

		JTable table = new JTable(data, columnNames);

		table.getColumn("transform").setCellRenderer(new CheckBoxRenderer());
		table.getColumn("transform").setCellEditor(new DefaultCellEditor(new JCheckBox()));

		table.getColumn("time invariant").setCellRenderer(new CheckBoxRenderer());
		table.getColumn("time invariant").setCellEditor(new DefaultCellEditor(new JCheckBox()));

		table.getColumn("remove predictor").setCellRenderer(new CheckBoxRenderer());
		table.getColumn("remove predictor").setCellEditor(new DefaultCellEditor(new JCheckBox()));

		JScrollPane scrollPane = new JScrollPane(table);
		scrollPane.setPreferredSize(new Dimension(450, 110));

		boxVert.add(scrollPane);

		// adds a new dummy variable to the covariates list
		updateButton.addActionListener(e -> {
			for (int i = covariateList.covariatesInput.get().size() - 1; i >= 0; i--) {
				if ((boolean) table.getModel().getValueAt(i, 3)) {
					covariateList.covariatesInput.get().remove(i);
				}else {
					covariateList.covariatesInput.get().get(i).logTransformInput.set(table.getModel().getValueAt(i, 1));
					covariateList.covariatesInput.get().get(i).timeInvariantInput.set(table.getModel().getValueAt(i, 2));
				}

			}
			refreshPanel();
		});

		add(boxVert);

	}

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

}
