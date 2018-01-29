/*
 * Copyright (C) 2015-2017 Jakub Galgonek   galgonek@uochb.cas.cz
 * Copyright (C) 2008-2009 Mark Rijnbeek    markr@ebi.ac.uk
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
package cz.iocb.sachem.shared;

import java.io.IOException;
import java.io.StringReader;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;



/**
 * Class for creating molecules
 */
public class MoleculeCreator
{
    private static final ThreadLocal<Aromaticity> aromaticity = new ThreadLocal<Aromaticity>()
    {
        @Override
        protected Aromaticity initialValue()
        {
            return new Aromaticity(ElectronDonation.cdkAllowingExocyclic(), Cycles.cdkAromaticSet());
        }
    };


    /**
     * Creates a molecule using an MDL string and an MDL reader .
     *
     * @param mdlString
     * @return
     * @throws CDKException
     * @throws IOException
     */
    public static IAtomContainer getMoleculeFromMolfile(String mol) throws CDKException, IOException
    {
        DefaultChemObjectReader mdlReader = null;

        if(mol.contains("M  V30 BEGIN CTAB"))
            mdlReader = new MDLV3000Reader();
        else
            mdlReader = new MDLV2000Reader();

        IAtomContainer readMolecule = new AtomContainer();
        mdlReader.setReader(new StringReader(mol));

        readMolecule = mdlReader.read(readMolecule);
        mdlReader.close();

        return readMolecule;
    }


    public static void configureMolecule(IAtomContainer molecule) throws CDKException
    {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        configureAromaticity(molecule);
    }


    public static void configureAromaticity(IAtomContainer molecule) throws CDKException
    {
        for(IBond bond : molecule.bonds())
            if(bond.isAromatic())
                return;

        aromaticity.get().apply(molecule);
    }
}
