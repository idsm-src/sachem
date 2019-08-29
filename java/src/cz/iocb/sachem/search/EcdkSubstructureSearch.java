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
package cz.iocb.sachem.search;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.BitSet;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeoutException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.sachem.fingerprint.ExtendedFingerprinter;
import cz.iocb.sachem.isomorphism.IsomorphismSort;
import cz.iocb.sachem.tautomers.CombinationCountException;
import cz.iocb.sachem.tautomers.InChIException;



public abstract class EcdkSubstructureSearch extends SubstructureSearch
{
    public static class EcdkQueryDataItem extends SubstructureSearch.QueryDataItem
    {
        //public short[] counts;
        public short[] fp;
    }


    private static final ThreadLocal<ExtendedFingerprinter> fingerPrinter = new ThreadLocal<ExtendedFingerprinter>()
    {
        @Override
        protected ExtendedFingerprinter initialValue()
        {
            return new ExtendedFingerprinter();
        }
    };


    public static QueryData<? extends EcdkQueryDataItem> getQueryData(byte[] queryArray, int type,
            boolean implicitHydrogens, boolean tautomers) throws CDKException, IOException, TimeoutException,
            CloneNotSupportedException, CombinationCountException, InChIException
    {
        String query = new String(queryArray, StandardCharsets.ISO_8859_1);

        List<IAtomContainer> queryMolecules;
        String message = null;

        try
        {
            queryMolecules = translateUserQuery(query, type, tautomers);
        }
        catch(CombinationCountException | InChIException | TimeoutException e)
        {
            queryMolecules = translateUserQuery(query, type, false);
            message = "cannot generate tautomers: " + e.getMessage();
        }
        catch(CDKException e)
        {
            queryMolecules = new LinkedList<IAtomContainer>();
            message = "cannot parse query '" + query + "'";
        }

        EcdkQueryDataItem[] items = new EcdkQueryDataItem[queryMolecules.size()];

        for(int idx = 0; idx < queryMolecules.size(); idx++)
        {
            IAtomContainer queryMolecule = queryMolecules.get(idx);

            IAtom[] sortedAtoms = IsomorphismSort.atomsByFrequency(queryMolecule);
            queryMolecule.setAtoms(sortedAtoms);

            BitSet fpbits = fingerPrinter.get().getFingerprint(queryMolecule);

            short[] fp = new short[fpbits.cardinality()];

            short p = 0;
            for(int i = fpbits.nextSetBit(0); i >= 0; i = fpbits.nextSetBit(i + 1))
                fp[p++] = (short) i;

            SachemMoleculeBuilder builder = new SachemMoleculeBuilder(queryMolecule);
            byte[] moleculeBytes = builder.asBytes(implicitHydrogens);

            //MoleculeCounts counts = new MoleculeCounts(queryMolecule, false);

            boolean[] restH = new boolean[queryMolecule.getAtomCount()];

            for(int i = 0; i < queryMolecule.getAtomCount(); i++)
            {
                IAtom atom = queryMolecule.getAtom(i);

                restH[i] = atom.getProperty(CDKConstants.REST_H) != null
                        && atom.getProperty(CDKConstants.REST_H).equals(true);
            }


            items[idx] = new EcdkQueryDataItem();

            /*
            data[idx].counts = new short[13];
            data[idx].counts[0] = counts.molSingleBondCount;
            data[idx].counts[1] = counts.molDoubleBondCount;
            data[idx].counts[2] = counts.molTripleBondCount;
            data[idx].counts[3] = counts.molAromaticBondCount;
            data[idx].counts[4] = counts.molSCount;
            data[idx].counts[5] = counts.molOCount;
            data[idx].counts[6] = counts.molNCount;
            data[idx].counts[7] = counts.molFCount;
            data[idx].counts[8] = counts.molClCount;
            data[idx].counts[9] = counts.molBrCount;
            data[idx].counts[10] = counts.molICount;
            data[idx].counts[11] = counts.molCCount;
            data[idx].counts[12] = counts.molPCount;
            */


            items[idx].fp = fp;
            items[idx].molecule = moleculeBytes;
            items[idx].restH = restH;
        }


        QueryData<EcdkQueryDataItem> data = new QueryData<EcdkQueryDataItem>();
        data.items = items;
        data.message = message;

        return data;
    }
}
