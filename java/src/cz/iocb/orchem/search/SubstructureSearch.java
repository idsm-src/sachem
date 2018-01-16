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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeoutException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.RGroupQueryReader;
import org.openscience.cdk.isomorphism.matchers.RGroupQuery;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import cz.iocb.orchem.convert.ConvertMolecule;
import cz.iocb.orchem.isomorphism.IsomorphismSort;
import cz.iocb.orchem.shared.MoleculeCreator;
import cz.iocb.orchem.tautomers.CombinationCountException;
import cz.iocb.orchem.tautomers.InchiTautomerGenerator;



public abstract class SubstructureSearch
{
    public static class QueryData
    {
        public byte[] molecule;
        public boolean[] restH;
    }


    private static final String QUERY_TYPE_MOL = "MOL";
    private static final String QUERY_TYPE_RGROUP = "RGROUP";
    private static final String QUERY_TYPE_SMILES = "SMILES";


    public static QueryData[] getQueryData(byte[] queryArray, String type, boolean implicitHydrogens, boolean tautomers)
            throws CDKException, IOException, TimeoutException, CloneNotSupportedException, CombinationCountException
    {
        String query = new String(queryArray, StandardCharsets.ISO_8859_1);
        List<IAtomContainer> queryMolecules = translateUserQuery(query, type, tautomers);
        QueryData[] data = new QueryData[queryMolecules.size()];

        for(int idx = 0; idx < queryMolecules.size(); idx++)
        {
            IAtomContainer queryMolecule = queryMolecules.get(idx);

            IAtom[] sortedAtoms = IsomorphismSort.atomsByFrequency(queryMolecule);
            queryMolecule.setAtoms(sortedAtoms);

            OrchemMoleculeBuilder builder = new OrchemMoleculeBuilder(queryMolecule);
            byte[] moleculeBytes = builder.asBytes(implicitHydrogens);

            boolean[] restH = new boolean[queryMolecule.getAtomCount()];

            for(int i = 0; i < queryMolecule.getAtomCount(); i++)
            {
                IAtom atom = queryMolecule.getAtom(i);

                restH[i] = atom.getProperty(CDKConstants.REST_H) != null
                        && atom.getProperty(CDKConstants.REST_H).equals(true);
            }

            data[idx] = new QueryData();
            data[idx].molecule = moleculeBytes;
            data[idx].restH = restH;
        }

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
     */
    protected static List<IAtomContainer> translateUserQuery(String query, String queryType, boolean tautomers)
            throws CDKException, IOException, TimeoutException, CloneNotSupportedException, CombinationCountException
    {
        List<IAtomContainer> userQueries = null;

        if(queryType.equals(QUERY_TYPE_RGROUP))
        {
            InputStream ins = new ByteArrayInputStream(query.getBytes());

            try(RGroupQueryReader reader = new RGroupQueryReader(ins))
            {
                RGroupQuery rGroupQuery = reader.read(new RGroupQuery(SilentChemObjectBuilder.getInstance()));
                userQueries = rGroupQuery.getAllConfigurations();
            }

            //TODO: add support for tautomers
        }
        else if(queryType.equals(QUERY_TYPE_MOL))
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
        else if(queryType.equals(QUERY_TYPE_SMILES))
        {
            if(!tautomers)
            {
                SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
                IAtomContainer molecule = sp.parseSmiles(query);

                userQueries = new ArrayList<IAtomContainer>();
                userQueries.add(molecule);
            }
            else
            {
                String molfile = ConvertMolecule.smilesToMolfile(query, false, false);
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
