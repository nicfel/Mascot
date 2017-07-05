/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package beast.mascot.beauti;

import beast.app.beauti.BeautiDoc;
import beast.app.beauti.GuessPatternDialog;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;

/**
 * BEAUti input editor for MultiTypeTree type traits.
 *
 * @author Nicola Felix Müller (nicola.felix.mueller@gmail.com)
 */
public class FilenameInputEditor extends InputEditor.Base {

    String filename;

    public FilenameInputEditor(BeautiDoc doc) {
        super(doc);

    }

    @Override
    public Class<?> type() {
        return String.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {

        JButton guessButton = new JButton("Search");
        guessButton.addActionListener((ActionEvent e) -> {
            refreshPanel();
        });

        JButton clearButton = new JButton("Type In");
        clearButton.addActionListener((ActionEvent e) -> {
            refreshPanel();
        });

        Box boxVert = Box.createVerticalBox();

        Box boxHoriz = Box.createHorizontalBox();
        boxHoriz.add(Box.createHorizontalGlue());
        boxHoriz.add(guessButton);
        boxHoriz.add(clearButton);
        boxVert.add(boxHoriz);

        add(boxVert);
    }
    
}
