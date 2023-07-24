package mascot.app.beauti;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beastfx.app.inputeditor.BeautiDoc;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Nicola Felix Müller <nicola.felix.muelelr@gmail.com>
 */
@Description("A phylogenetic tree with a trait class (somehow needed for Beauti)")
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

            //typeList = Lists.newArrayList(typeSet);
            typeList = new ArrayList<>();
            typeList.addAll(typeSet);

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
