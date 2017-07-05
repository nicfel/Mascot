/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import com.google.common.collect.Lists;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("A multi-type phylogenetic tree.")
@Citation("Timothy G. Vaughan, Denise Kuhnert, Alex Popinga, David Welch and \n"
        + "Alexei J. Drummond, 'Efficient Bayesian inference under the \n"
        + "structured coalescent', Bioinformatics 30:2272, 2014.")
public class TreeWithTrait extends Tree implements StateNodeInitialiser {

    /*
     * Inputs:
     */
    public Input<String> typeLabelInput = new Input<>(
        "typeLabel",
        "Label for type traits (default 'type')", "type");

    public Input<TraitSet> typeTraitInput = new Input<>(
        "typeTrait", "Type trait set.  Used only by BEAUti.");

    public Input<String> typeTraitValuesInput = new Input<>(
            "typeTraitValues",
            "Comma-delimited list of types to be included even when absent " +
                    "from the sampled taxa.");

    /*
     * Non-input fields:
     */
    protected String typeLabel;
    protected TraitSet typeTraitSet;
    
    protected List <String> typeList;

    
    @Override
    public void initAndValidate() {
    	super.initAndValidate();
        typeLabel = typeLabelInput.get();
        
        processTraits(m_traitList.get());

        // Ensure tree is compatible with traits.
        if (hasDateTrait())
            adjustTreeNodeHeights(root);
    }

    @Override
    protected void processTraits(List<TraitSet> traitList) {
        super.processTraits(traitList);
        
        // Record trait set associated with leaf types.
        for (TraitSet traitSet : traitList) {
            if (traitSet.getTraitName().equals(typeLabel)) {
                typeTraitSet = traitSet;
                break;
            }
        }

        // Use explicitly-identified type trait set if available.
        // Seems dumb, but needed for BEAUti as ListInputEditors
        // muck things up...
        if (typeTraitInput.get() != null)
            typeTraitSet = typeTraitInput.get();

        // Construct type list.
        if (typeTraitSet == null) {
            if (getTaxonset() != null) {
                TraitSet dummyTraitSet = new TraitSet();

                StringBuilder sb = new StringBuilder();
                for (int i=0; i<getTaxonset().getTaxonCount(); i++) {
                    if (i>0)
                        sb.append(",\n");
                    sb.append(getTaxonset().getTaxonId(i)).append("=NOT_SET");
                }
                try {
                    dummyTraitSet.initByName(
                        "traitname", "type",
                        "taxa", getTaxonset(),
                        "value", sb.toString());
                    dummyTraitSet.setID("typeTraitSet.t:"
                        + BeautiDoc.parsePartition(getID()));
                    setTypeTrait(dummyTraitSet);
                } catch (Exception ex) {
                    System.out.println("Error setting default type trait.");
                }
            }
        }

        if (typeTraitSet != null) {

            Set<String> typeSet = new HashSet<>();

            int nTaxa = typeTraitSet.taxaInput.get().asStringList().size();
            for (int i = 0; i < nTaxa; i++)
                typeSet.add(typeTraitSet.getStringValue(i));

            // Include any addittional trait values in type list
            if (typeTraitValuesInput.get() != null) {
                for (String typeName : typeTraitValuesInput.get().split(","))
                    typeSet.add(typeName);
            }

            typeList = Lists.newArrayList(typeSet);
            Collections.sort(typeList);

            System.out.println("Type trait with the following types detected:");
            for (int i = 0; i < typeList.size(); i++)
                System.out.println(typeList.get(i) + " (" + i + ")");

        }
    }
    
    /**
     * @return TraitSet with same name as typeLabel.
     */
    public TraitSet getTypeTrait() {
        if (!traitsProcessed)
            processTraits(m_traitList.get());
        
        return typeTraitSet;
    }
    
    /**
     * @return true if TraitSet with same name as typeLabel exists.
     */
    public boolean hasTypeTrait() {
        return getTypeTrait() != null;
    }

    /**
     * Specifically set the type trait set for this tree. A null value simply
     * removes the existing trait set.
     *
     * @param traitSet
     */
    public void setTypeTrait(TraitSet traitSet) {
        if (hasTypeTrait()) {
            m_traitList.get().remove(typeTraitSet);
        }

        if (traitSet != null) {
            //m_traitList.setValue(traitSet, this);
            typeTraitInput.setValue(traitSet, this);
        }

        typeTraitSet = traitSet;
    }
    
    /**
     * Retrieve the list of unique types identified by the type trait.
     * @return List of unique type trait value strings.
     */
    public List<String> getTypeList() {
        if (!traitsProcessed)
            processTraits(m_traitList.get());
        
        return typeList;
    }

    /**
     * @param type
     * @return string name of given type
     */
    public String getTypeString(int type) {
        if (!traitsProcessed)
            processTraits(m_traitList.get());

        return typeList.get(type);
    }

    /**
     * @param typeString
     * @return integer type corresponding to given type string
     */
    public int getTypeFromString(String typeString) {
        if (!traitsProcessed)
            processTraits(m_traitList.get());

        return typeList.indexOf(typeString);
    }

    /**
     * @return type label to be used in logging.
     */
    public String getTypeLabel() {
        return typeLabel;
    }


    @Override
    public void initStateNodes() { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodeList) {
        stateNodeList.add(this);
    }    



}
