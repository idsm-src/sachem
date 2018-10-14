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
import java.util.LinkedList;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;



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


        /* fix deuterium and tritium pseudo atoms */
        LinkedList<IAtom> pseudoHydrogens = new LinkedList<IAtom>();
        LinkedList<IAtom> newHydrogens = new LinkedList<IAtom>();

        for(IAtom atom : readMolecule.atoms())
            if(atom instanceof IPseudoAtom && (atom.getSymbol().equals("D") || atom.getSymbol().equals("T")))
                pseudoHydrogens.add(atom);

        IChemObjectBuilder builder = readMolecule.getBuilder();

        for(IAtom atom : pseudoHydrogens)
        {
            IIsotope isotope = builder.newInstance(IIsotope.class, "H", atom.getSymbol().equals("D") ? 2 : 3);
            IAtom hydrogen = builder.newInstance(IAtom.class, isotope);

            newHydrogens.add(hydrogen);
            readMolecule.addAtom(hydrogen);

            for(IBond bond : readMolecule.bonds())
            {
                if(bond.contains(atom))
                {
                    IAtom other = bond.getOther(atom);
                    IAtom[] atoms = { other, hydrogen };
                    bond.setAtoms(atoms);
                }
            }

            readMolecule.removeAtom(atom);
        }


        CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(readMolecule.getBuilder());

        for(IAtom atom : newHydrogens)
        {
            IAtomType type = matcher.findMatchingAtomType(readMolecule, atom);
            AtomTypeManipulator.configure(atom, type);
            CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(readMolecule.getBuilder());
            adder.addImplicitHydrogens(readMolecule, atom);
        }

        for(IAtom atom : readMolecule.atoms())
        {
            if(atom.getImplicitHydrogenCount() == null)
            {
                IAtomType type = matcher.findMatchingAtomType(readMolecule, atom);
                AtomTypeManipulator.configure(atom, type);
                CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(readMolecule.getBuilder());
                adder.addImplicitHydrogens(readMolecule, atom);
            }
        }

        return readMolecule;
    }


    public static IAtomContainer getMoleculeFromSmiles(String smiles) throws CDKException
    {
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        sp.kekulise(false);
        IAtomContainer molecule = sp.parseSmiles(smiles);

        return molecule;
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
