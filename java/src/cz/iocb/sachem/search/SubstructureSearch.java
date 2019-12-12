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
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeoutException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.sachem.molecule.AromaticityMode;
import cz.iocb.sachem.molecule.BinaryMoleculeBuilder;
import cz.iocb.sachem.molecule.MoleculeCreator;
import cz.iocb.sachem.molecule.TautomerMode;
import cz.iocb.sachem.tautomers.CombinationCountException;
import cz.iocb.sachem.tautomers.InChIException;



public abstract class SubstructureSearch
{
    public static class QueryData<ItemType extends QueryDataItem>
    {
        ItemType[] items;
        String message;
    }


    public static class QueryDataItem
    {
        public byte[] molecule;
        public boolean[] restH;
    }


    public static QueryData<? extends QueryDataItem> getQueryData(byte[] queryArray, int type,
            boolean implicitHydrogens, boolean tautomers) throws CDKException, IOException, TimeoutException,
            CloneNotSupportedException, CombinationCountException, InChIException
    {
        String query = new String(queryArray, StandardCharsets.ISO_8859_1);

        List<IAtomContainer> queryMolecules = null;
        String message = null;

        try
        {
            queryMolecules = MoleculeCreator.translateQuery(query, AromaticityMode.PRESERVE,
                    tautomers ? TautomerMode.INCHI : TautomerMode.IGNORE);
        }
        catch(CombinationCountException | InChIException | TimeoutException e)
        {
            queryMolecules = MoleculeCreator.translateQuery(query, AromaticityMode.PRESERVE, TautomerMode.IGNORE);
            message = "cannot generate tautomers: " + e.getMessage();
        }
        catch(CDKException e)
        {
            queryMolecules = new LinkedList<IAtomContainer>();
            message = "cannot parse query '" + query + "'";
        }

        QueryDataItem[] items = new QueryDataItem[queryMolecules.size()];

        for(int idx = 0; idx < queryMolecules.size(); idx++)
        {
            IAtomContainer queryMolecule = queryMolecules.get(idx);

            BinaryMoleculeBuilder builder = new BinaryMoleculeBuilder(queryMolecule);
            byte[] moleculeBytes = builder.asBytes(implicitHydrogens);

            boolean[] restH = new boolean[queryMolecule.getAtomCount()];

            for(int i = 0; i < queryMolecule.getAtomCount(); i++)
            {
                IAtom atom = queryMolecule.getAtom(i);

                restH[i] = atom.getProperty(CDKConstants.REST_H) != null
                        && atom.getProperty(CDKConstants.REST_H).equals(true);
            }

            items[idx] = new QueryDataItem();
            items[idx].molecule = moleculeBytes;
            items[idx].restH = restH;
        }


        QueryData<QueryDataItem> data = new QueryData<QueryDataItem>();
        data.items = items;
        data.message = message;

        return data;
    }
}
