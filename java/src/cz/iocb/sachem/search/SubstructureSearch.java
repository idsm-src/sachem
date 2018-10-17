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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeoutException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.RGroupQueryReader;
import org.openscience.cdk.isomorphism.matchers.RGroupQuery;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import cz.iocb.sachem.isomorphism.IsomorphismSort;
import cz.iocb.sachem.shared.MoleculeCreator;
import cz.iocb.sachem.tautomers.CombinationCountException;
import cz.iocb.sachem.tautomers.InChIException;
import cz.iocb.sachem.tautomers.InchiTautomerGenerator;



public abstract class SubstructureSearch
{
    public static class QueryData
    {
        QueryDataItem[] items;
        String message;
    }


    public static class QueryDataItem
    {
        public byte[] molecule;
        public boolean[] restH;
    }


    public static QueryData getQueryData(byte[] queryArray, int type, boolean implicitHydrogens, boolean tautomers)
            throws CDKException, IOException, TimeoutException, CloneNotSupportedException, CombinationCountException,
            InChIException
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

        QueryDataItem[] items = new QueryDataItem[queryMolecules.size()];

        for(int idx = 0; idx < queryMolecules.size(); idx++)
        {
            IAtomContainer queryMolecule = queryMolecules.get(idx);

            IAtom[] sortedAtoms = IsomorphismSort.atomsByFrequency(queryMolecule);
            queryMolecule.setAtoms(sortedAtoms);

            SachemMoleculeBuilder builder = new SachemMoleculeBuilder(queryMolecule);
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


        QueryData data = new QueryData();
        data.items = items;
        data.message = message;

        return data;
    }


    /**
     * Translates a user query (such as for example a Smiles string "O=S=O") into a CDK atomcontainer.
     *
     * However, if the query is an RGFile (a special case of a MOL file), there can be any number of queries generated
     * (valid substitutes of the R-Group)
     *
     * @param userQuery mdl or smiles clob
     * @param queryType SMILES or MOL
     * @param tautomers Y/N indicating to query for tautomers
     * @return user query represented as CDK atom container
     * @throws CDKException
     * @throws IOException
     * @throws CombinationCountException
     * @throws CloneNotSupportedException
     * @throws TimeoutException
     * @throws InChIException
     */
    protected static List<IAtomContainer> translateUserQuery(String query, int queryTypeIdx, boolean tautomers)
            throws CDKException, IOException, TimeoutException, CloneNotSupportedException, CombinationCountException,
            InChIException
    {
        QueryFormat queryType = QueryFormat.values()[queryTypeIdx];
        List<IAtomContainer> userQueries = null;

        if(queryType == QueryFormat.UNSPECIFIED)
            queryType = QueryFormat.detect(query);


        if(queryType == QueryFormat.RGROUP)
        {
            InputStream ins = new ByteArrayInputStream(query.getBytes());

            try(RGroupQueryReader reader = new RGroupQueryReader(ins))
            {
                RGroupQuery rGroupQuery = reader.read(new RGroupQuery(SilentChemObjectBuilder.getInstance()));
                userQueries = rGroupQuery.getAllConfigurations();
            }

            if(tautomers)
            {
                List<IAtomContainer> queries = new ArrayList<IAtomContainer>();

                for(IAtomContainer molecule : userQueries)
                {
                    String molfile = MoleculeCreator.getMolfileFromMolecule(molecule);
                    InchiTautomerGenerator tautomerGenerator = new InchiTautomerGenerator();
                    queries.addAll(tautomerGenerator.getTautomers(molfile));
                }

                userQueries = queries;
            }
        }
        else if(queryType == QueryFormat.MOLFILE)
        {
            if(!tautomers)
            {
                IAtomContainer readMolecule = MoleculeCreator.getMoleculeFromMolfile(query);

                userQueries = new ArrayList<IAtomContainer>();
                userQueries.add(readMolecule);
            }
            else
            {
                InchiTautomerGenerator tautomerGenerator = new InchiTautomerGenerator();
                userQueries = tautomerGenerator.getTautomers(query);
            }
        }
        else if(queryType == QueryFormat.SMILES)
        {
            IAtomContainer molecule = MoleculeCreator.getMoleculeFromSmiles(query);

            if(!tautomers)
            {
                userQueries = new ArrayList<IAtomContainer>();
                userQueries.add(molecule);
            }
            else
            {
                String molfile = MoleculeCreator.getMolfileFromMolecule(molecule);
                InchiTautomerGenerator tautomerGenerator = new InchiTautomerGenerator();
                userQueries = tautomerGenerator.getTautomers(molfile);
            }
        }
        else
        {
            userQueries = new ArrayList<IAtomContainer>(0);
        }


        if(userQueries != null)
            for(IAtomContainer molecule : userQueries)
                MoleculeCreator.configureMolecule(molecule);


        return userQueries;
    }
}
