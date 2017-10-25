/*
 * Copyright (C) 2015-2017 Jakub Galgonek   galgonek@uochb.cas.cz
 * Copyright (C) 2008-2010 Mark Rijnbeek    markr@ebi.ac.uk
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
package cz.iocb.orchem.search;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.BitSet;
import java.util.List;
import java.util.concurrent.TimeoutException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.orchem.fingerprint.OrchemFingerprinter;
import cz.iocb.orchem.isomorphism.IsomorphismSort;
import cz.iocb.orchem.shared.MoleculeCounts;
import cz.iocb.orchem.tautomers.CombinationCountException;



public abstract class OrchemSubstructureSearch extends SubstructureSearch
{
    public static class OrchemQueryData extends SubstructureSearch.QueryData
    {
        public short[] counts;
        public short[] fp;
    }


    private static final ThreadLocal<OrchemFingerprinter> fingerPrinter = new ThreadLocal<OrchemFingerprinter>()
    {
        @Override
        protected OrchemFingerprinter initialValue()
        {
            return new OrchemFingerprinter();
        }
    };


    public static OrchemQueryData[] getQueryData(byte[] queryArray, String type, boolean tautomers)
            throws CDKException, IOException, TimeoutException, CloneNotSupportedException, CombinationCountException
    {
        String query = new String(queryArray, StandardCharsets.ISO_8859_1);
        List<IAtomContainer> queryMolecules = translateUserQuery(query, type, tautomers);
        OrchemQueryData[] data = new OrchemQueryData[queryMolecules.size()];

        for(int idx = 0; idx < queryMolecules.size(); idx++)
        {
            IAtomContainer queryMolecule = queryMolecules.get(idx);

            IAtom[] sortedAtoms = IsomorphismSort.atomsByFrequency(queryMolecule);
            queryMolecule.setAtoms(sortedAtoms);

            BitSet fpbits = fingerPrinter.get().getBitFingerprint(queryMolecule).asBitSet();

            short[] fp = new short[fpbits.cardinality() - 1];

            short p = 0;
            for(int i = fpbits.nextSetBit(1); i >= 0; i = fpbits.nextSetBit(i + 1))
                fp[p++] = (short) (i - 1);

            OrchemMoleculeBuilder builder = new OrchemMoleculeBuilder(queryMolecule);
            byte[] atomBytes = builder.atomsAsBytes();
            byte[] bondBytes = builder.bondsAsBytes();

            MoleculeCounts counts = new MoleculeCounts(queryMolecule);

            boolean[] restH = new boolean[queryMolecule.getAtomCount()];

            for(int i = 0; i < queryMolecule.getAtomCount(); i++)
            {
                IAtom atom = queryMolecule.getAtom(i);

                restH[i] = atom.getProperty(CDKConstants.REST_H) != null
                        && atom.getProperty(CDKConstants.REST_H).equals(true);
            }


            data[idx] = new OrchemQueryData();

            data[idx].counts = new short[10];
            data[idx].counts[0] = counts.molTripleBondCount;
            data[idx].counts[1] = counts.molSCount;
            data[idx].counts[2] = counts.molOCount;
            data[idx].counts[3] = counts.molNCount;
            data[idx].counts[4] = counts.molFCount;
            data[idx].counts[5] = counts.molClCount;
            data[idx].counts[6] = counts.molBrCount;
            data[idx].counts[7] = counts.molICount;
            data[idx].counts[8] = counts.molCCount;
            data[idx].counts[9] = counts.molPCount;

            data[idx].fp = fp;
            data[idx].atoms = atomBytes;
            data[idx].bonds = bondBytes;
            data[idx].restH = restH;
        }

        return data;
    }
}
