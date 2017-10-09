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
package cz.iocb.orchem.shared;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry.Conformation;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.ExtendedTetrahedral;
import org.openscience.cdk.stereo.TetrahedralChirality;
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
    public static AtomContainer getMoleculeFromMolfile(String mol) throws CDKException, IOException
    {
        DefaultChemObjectReader mdlReader = null;

        if(mol.contains("M  V30 BEGIN CTAB"))
            mdlReader = new MDLV3000Reader();
        else
            mdlReader = new MDLV2000Reader();

        AtomContainer readMolecule = new AtomContainer();
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
        aromaticity.get().apply(molecule);
    }


    /**
     * HACK HACK - should be taken from CDK but need a fix first for the builder, because we can't get a
     * DefaultChemObjectBuilder from a query atom container. Produces an AtomContainer without explicit Hs but with H
     * count from one with Hs. The new molecule is a deep copy.
     *
     * @param atomContainer The AtomContainer from which to remove the hydrogens
     * @return The molecule without Hs.
     * @cdk.keyword hydrogens, removal
     */
    public static IAtomContainer removeHydrogens(IAtomContainer atomContainer)
    {
        Map<IAtom, IAtom> map = new HashMap<IAtom, IAtom>(); // maps original atoms to clones
        List<IAtom> remove = new ArrayList<IAtom>(); // lists removed hydrogens

        IChemObjectBuilder bob = SilentChemObjectBuilder.getInstance();
        IAtomContainer mol = bob.newInstance(IAtomContainer.class);

        int count = atomContainer.getAtomCount();
        for(int i = 0; i < count; i++)
        {
            IAtom atom = atomContainer.getAtom(i);

            if(!atom.getSymbol().equals("H"))
            {
                IAtom clonedAtom = null;
                try
                {
                    clonedAtom = atom.clone();
                }
                catch (CloneNotSupportedException e)
                {
                    e.printStackTrace();
                }
                mol.addAtom(clonedAtom);
                map.put(atom, clonedAtom);
            }
            else
            {
                remove.add(atom);
            }
        }

        count = atomContainer.getBondCount();
        for(int i = 0; i < count; i++)
        {
            final IBond bond = atomContainer.getBond(i);
            boolean removedBond = false;
            final int length = bond.getAtomCount();
            for(int k = 0; k < length; k++)
            {
                if(remove.contains(bond.getAtom(k)))
                {
                    removedBond = true;
                    break;
                }
            }
            if(!removedBond)
            {
                IBond clone = null;
                try
                {
                    clone = atomContainer.getBond(i).clone();
                }
                catch (CloneNotSupportedException e)
                {
                    e.printStackTrace();
                }
                assert clone != null;
                clone.setAtoms(new IAtom[] { map.get(bond.getAtom(0)), map.get(bond.getAtom(1)) });
                mol.addBond(clone);
            }
        }

        for(IAtom aRemove : remove)
        {
            for(IAtom iAtom : atomContainer.getConnectedAtomsList(aRemove))
            {
                final IAtom neighb = map.get(iAtom);

                if(neighb == null)
                    continue; // since for the case of H2, neight H has a heavy atom neighbor

                neighb.setImplicitHydrogenCount(
                        (neighb.getImplicitHydrogenCount() == null ? 0 : neighb.getImplicitHydrogenCount()) + 1);
            }
        }

        for(IAtom atom : mol.atoms())
        {
            if(atom.getImplicitHydrogenCount() == null)
                atom.setImplicitHydrogenCount(0);
        }

        mol.setProperties(atomContainer.getProperties());
        mol.setFlags(atomContainer.getFlags());


        for(IStereoElement element : atomContainer.stereoElements())
        {
            if(element instanceof TetrahedralChirality)
            {
                TetrahedralChirality org = (TetrahedralChirality) element;
                IAtom[] orgLigands = org.getLigands();

                IAtom chiralAtom = map.get(org.getChiralAtom());
                IAtom[] ligandAtoms = new IAtom[orgLigands.length];

                for(int i = 0; i < orgLigands.length; i++)
                {
                    ligandAtoms[i] = map.get(orgLigands[i]);

                    if(ligandAtoms[i] == null)
                        ligandAtoms[i] = chiralAtom;
                }

                mol.addStereoElement(new TetrahedralChirality(chiralAtom, ligandAtoms, org.getStereo()));
            }
            else if(element instanceof ExtendedTetrahedral)
            {
                ExtendedTetrahedral org = (ExtendedTetrahedral) element;
                IAtom[] orgPeripherals = org.peripherals();

                IAtom focusAtom = map.get(org.focus());
                IAtom[] peripheralsAtoms = new IAtom[orgPeripherals.length];

                for(int i = 0; i < orgPeripherals.length; i++)
                {
                    peripheralsAtoms[i] = map.get(orgPeripherals[i]);

                    if(peripheralsAtoms[i] == null)
                        peripheralsAtoms[i] = map.get(org.findTerminalAtoms(atomContainer)[i < 2 ? 0 : 1]);
                }

                mol.addStereoElement(new ExtendedTetrahedral(focusAtom, peripheralsAtoms, org.winding()));
            }
            else if(element instanceof DoubleBondStereochemistry)
            {
                DoubleBondStereochemistry org = (DoubleBondStereochemistry) element;

                IAtom stereoAtom0 = map.get(org.getStereoBond().getAtom(0));
                IAtom stereoAtom1 = map.get(org.getStereoBond().getAtom(1));
                IBond stereoBond = mol.getBond(stereoAtom0, stereoAtom1);

                IBond[] orgBonds = org.getBonds();

                IBond[] bonds = new IBond[orgBonds.length];
                Conformation stereo = org.getStereo();

                for(int i = 0; i < orgBonds.length; i++)
                {
                    IAtom a0 = map.get(orgBonds[i].getAtom(0));
                    IAtom a1 = map.get(orgBonds[i].getAtom(1));

                    if(a0 == null)
                    {
                        List<IAtom> list = mol.getConnectedAtomsList(a1);
                        list.remove(stereoBond.getOther(a1));
                        a0 = list.get(0);

                        if(stereo == Conformation.OPPOSITE)
                            stereo = Conformation.TOGETHER;
                        else
                            stereo = Conformation.OPPOSITE;
                    }
                    else if(a1 == null)
                    {
                        List<IAtom> list = mol.getConnectedAtomsList(a0);
                        list.remove(stereoBond.getOther(a0));
                        a1 = list.get(0);

                        if(stereo == Conformation.OPPOSITE)
                            stereo = Conformation.TOGETHER;
                        else
                            stereo = Conformation.OPPOSITE;
                    }

                    bonds[i] = mol.getBond(a0, a1);
                }

                mol.addStereoElement(new DoubleBondStereochemistry(stereoBond, bonds, stereo));
            }
        }

        return mol;
    }
}
