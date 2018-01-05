/*
 * Copyright (C) 2015-2017 Jakub Galgonek   galgonek@uochb.cas.cz
 * Copyright (C) 2009-2009 Mark Rijnbeek    markr@ebi.ac.uk
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work, which includes - but is not limited to -
 * adding the above copyright notice to the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 */
package cz.iocb.orchem.isomorphism;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;



/**
 * Facilitates sort before doing an isomorphic comparison.
 */
public class IsomorphismSort
{
    /**
     * Returns a sorted array of {@link IAtom} for a given input {@link IAtomContainer}. The sort is based on:
     * <ul>
     * <li>the frequency of atom symbols</li>
     * <li>the connectivity of the atoms</li>
     * </ul>
     *
     * This type of sort benefits algorithms like VF2, sorting the query container before a subgraph match. The
     * algorithm will then start with atoms that are quite likely hardest to match.
     *
     * Example usage:
     *
     * <pre>
     * IAtom[] sortedAtoms = cz.iocb.orchem.IsomorphismSort.atomsByFrequency(queryMolecule);
     * queryMolecule.setAtoms(sortedAtoms);
     * </pre>
     *
     * @param iac input {@link IAtomContainer}
     * @return a sorted array of atoms for that {@link IAtomContainer}
     */
    public static IAtom[] atomsByFrequency(IAtomContainer iac)
    {
        // Create a map with (key,value) being (atom symbol, overall count)
        Map<String, Integer> elementCounts = new TreeMap<String, Integer>();

        for(IAtom atom : iac.atoms())
        {
            Integer count = elementCounts.get(atom.getSymbol());

            if(count == null)
                elementCounts.put(atom.getSymbol(), new Integer(1));
            else
                elementCounts.put(atom.getSymbol(), ++count);
        }

        // Create a map with (key,value) being (IAtom, number of bond IAtom occurs in)
        Map<IAtom, Integer> bondParticipationCount = new HashMap<IAtom, Integer>();

        for(IBond bond : iac.bonds())
        {
            if(bond.getAtom(0).getSymbol().equals("H") || bond.getAtom(1).getSymbol().equals("H"))
                continue;

            for(IAtom atomInBond : bond.atoms())
            {
                if(!bondParticipationCount.containsKey(atomInBond))
                    bondParticipationCount.put(atomInBond, 1);
                else
                    bondParticipationCount.put(atomInBond, bondParticipationCount.get(atomInBond) + 1);
            }
        }

        // Mind the atoms not in any bond
        for(IAtom atom : iac.atoms())
        {
            if(!bondParticipationCount.containsKey(atom))
                bondParticipationCount.put(atom, 0);
        }

        // We now have to maps that will be used to sort the incoming atom container
        List<AtomForIsomorphismSort> atomList = new ArrayList<AtomForIsomorphismSort>();

        for(IAtom atom : iac.atoms())
        {
            AtomForIsomorphismSort afis = new AtomForIsomorphismSort(atom, elementCounts.get(atom.getSymbol()),
                    bondParticipationCount.get(atom));
            atomList.add(afis);
        }

        Collections.sort(atomList, new AtomForIsomorphismSortComparator());

        // Create an output atom array based on the sorted list
        IAtom[] iAtomArray = new IAtom[iac.getAtomCount()];
        int iacSortedIdx = 0;

        for(AtomForIsomorphismSort afis : atomList)
        {
            iAtomArray[iacSortedIdx] = afis.iatom;
            iacSortedIdx++;
        }

        return iAtomArray;
    }


    /**
     * Beans helper class for sorting. It holds the values to sort, and has comparator class
     * {@link AtomForIsomorphismSortComparator} to determine the right sort order for these.
     */
    private static class AtomForIsomorphismSort
    {
        IAtom iatom;
        Integer overallElementCount;
        Integer atomBondParticipationCount;

        public AtomForIsomorphismSort(IAtom _iatom, int _overallElementCount, int _atomBondParticipationCount)
        {
            iatom = _iatom;
            overallElementCount = _overallElementCount;
            atomBondParticipationCount = _atomBondParticipationCount;
        }
    }


    /**
     * Comparator for the sort. The sort is based primary on the elements with lowest count first, and secondary on the
     * atoms that are most highly connected first.
     */
    private static class AtomForIsomorphismSortComparator implements Comparator<AtomForIsomorphismSort>
    {
        @Override
        public int compare(AtomForIsomorphismSort e1, AtomForIsomorphismSort e2)
        {
            if(e1.iatom.getSymbol().equals("H") && !e2.iatom.getSymbol().equals("H"))
                return 1;

            if(!e1.iatom.getSymbol().equals("H") && e2.iatom.getSymbol().equals("H"))
                return -1;

            if(e1.overallElementCount.compareTo(e2.overallElementCount) == 0)
                return e2.atomBondParticipationCount.compareTo(e1.atomBondParticipationCount);
            else
                return e1.overallElementCount.compareTo(e2.overallElementCount);
        }
    }
}
