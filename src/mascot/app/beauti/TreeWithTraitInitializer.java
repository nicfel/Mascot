/*
 * Copyright (C) 2012 Tim Vaughan
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
package mascot.app.beauti;

import beast.base.core.Description;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;

import java.util.List;

/**
 * A multi-type tree generated randomly from leaf types and a migration matrix
 * with fixed population sizes.
 *
 * @author Tim Vaughan
 */
@Description("A multi-type tree generated randomly from leaf types and"
+ "a migration matrix with fixed population sizes.")
public class TreeWithTraitInitializer extends TreeWithTrait implements StateNodeInitialiser {
   

    
    public TreeWithTraitInitializer() { }

    @Override
    public void initAndValidate() {
    }

    
    @Override
    public void initStateNodes() { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodeList) {
        stateNodeList.add(this);
    }    

}
