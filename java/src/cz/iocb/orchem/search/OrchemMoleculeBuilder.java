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
package cz.iocb.orchem.search;

import java.util.List;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry.Conformation;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.ExtendedTetrahedral;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.stereo.Stereocenters.Type;
import org.openscience.cdk.stereo.TetrahedralChirality;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;



public class OrchemMoleculeBuilder
{
    private static final int atomBlockSize = 3;
    private static final int bondBlockSize = 4;
    private static final int hAtomBlockSize = 1;
    private static final int hBondBlockSize = 3;

    private static int[] validReorder = { 0x1234, 0x1423, 0x1342, 0x2314, 0x2431, 0x2143, 0x3124, 0x3412, 0x3241,
            0x4213, 0x4321, 0x4132 };

    private IAtomContainer molecule;
    private Stereocenters centers;


    public OrchemMoleculeBuilder(IAtomContainer molecule)
    {
        this.molecule = molecule;
        centers = Stereocenters.of(molecule);
    }


    public byte[] atomsAsBytes()
    {
        int count = 0;

        while(count < molecule.getAtomCount() && !molecule.getAtom(count).getSymbol().equals("H"))
            count++;

        for(int i = count; i < molecule.getAtomCount(); i++)
            assert molecule.getAtom(i).getSymbol().equals("H");


        byte[] array = new byte[count * atomBlockSize];
        int offset = 0;

        for(IAtom atom : molecule.atoms())
        {
            String symbol = atom.getSymbol();

            if(atom.getSymbol().equals("H"))
                continue;


            if(atom instanceof IPseudoAtom)
            {
                String label = ((IPseudoAtom) atom).getLabel();
                assert label.length() == 0 || label.length() > 1 || label.charAt(0) >= 0 && label.charAt(0) < 128;

                if(label.length() == 1)
                    array[offset + 0] = (byte) -label.charAt(0);
                else if(label.matches("R[1-9#]?"))
                    array[offset + 0] = -'R';
                else
                    array[offset + 0] = -'?';
            }
            else
            {
                int num = PeriodicTable.getAtomicNumber(symbol);
                assert num > 0 && num < 128;

                array[offset + 0] = (byte) num;
            }


            Integer implicitH = atom.getImplicitHydrogenCount();
            int explicitH = AtomContainerManipulator.countExplicitHydrogens(molecule, atom);

            if(implicitH == null)
                implicitH = 0;

            assert implicitH <= 31;
            assert explicitH <= 31;

            array[offset + 1] = (byte) (explicitH << 4 | implicitH);
            array[offset + 2] = (byte) (atom.getFormalCharge() << 2);

            offset += atomBlockSize;
        }


        for(IStereoElement element : molecule.stereoElements())
        {
            if(element instanceof TetrahedralChirality)
            {
                TetrahedralChirality chirality = (TetrahedralChirality) element;
                int flag = getChiralityValue(chirality.getChiralAtom(), chirality.getLigands(), chirality.getStereo());

                int index = molecule.indexOf(chirality.getChiralAtom());
                array[index * atomBlockSize + 2] |= flag;
            }
            else if(element instanceof ExtendedTetrahedral)
            {
                ExtendedTetrahedral chirality = (ExtendedTetrahedral) element;
                int flag = getChiralityValue(chirality.focus(), chirality.peripherals(), chirality.winding());

                int index = molecule.indexOf(chirality.focus());
                array[index * atomBlockSize + 2] |= flag;
            }
        }


        for(int index = 0; index < count; index++)
            if((array[index * atomBlockSize + 2] & 0x03) == 0
                    && (isTetrahedralChirality(index) || isExtendedTetrahedral(index)))
                array[index * atomBlockSize + 2] |= 0x03;


        return array;
    }


    public byte[] hAtomsAsBytes()
    {
        int nonH = 0;

        while(nonH < molecule.getAtomCount() && !molecule.getAtom(nonH).getSymbol().equals("H"))
            nonH++;

        for(int i = nonH; i < molecule.getAtomCount(); i++)
            assert molecule.getAtom(i).getSymbol().equals("H");

        int count = molecule.getAtomCount() - nonH;


        byte[] array = new byte[count * hAtomBlockSize];
        int offset = 0;

        for(IAtom atom : molecule.atoms())
        {
            if(!atom.getSymbol().equals("H"))
                continue;

            array[offset + 0] = (byte) (int) atom.getFormalCharge();
            offset += hAtomBlockSize;
        }

        return array;
    }


    public byte[] bondsAsBytes()
    {
        DoubleBondStereochemistry[] stereochemistry = new DoubleBondStereochemistry[molecule.getBondCount()];

        for(IStereoElement element : molecule.stereoElements())
        {
            if(element instanceof DoubleBondStereochemistry)
            {
                DoubleBondStereochemistry e = (DoubleBondStereochemistry) element;
                stereochemistry[molecule.indexOf(e.getStereoBond())] = e;
            }
        }


        int count = 0;

        for(IBond bond : molecule.bonds())
            if(!bond.getAtom(0).getSymbol().equals("H") && !bond.getAtom(1).getSymbol().equals("H"))
                count++;


        byte[] array = new byte[count * bondBlockSize];
        int offset = 0;

        for(int idx = 0; idx < molecule.getBondCount(); idx++)
        {
            IBond bond = molecule.getBond(idx);
            assert bond.getAtomCount() == 2;

            IAtom a1 = bond.getAtom(0);
            IAtom a2 = bond.getAtom(1);

            if(a1.getSymbol().equals("H") || a2.getSymbol().equals("H"))
                continue;


            int a1idx = molecule.indexOf(a1);
            int a2idx = molecule.indexOf(a2);
            assert a1idx != a2idx;

            array[offset + 0] = (byte) (a1idx % 256);
            array[offset + 1] = (byte) (a1idx / 256 << 4 | a2idx / 256);
            array[offset + 2] = (byte) (a2idx % 256);


            assert bond.getOrder().ordinal() < 32;
            array[offset + 3] = (byte) (bond.getOrder().ordinal() << 2);
            array[offset + 3] |= (byte) (bond.getFlag(CDKConstants.ISAROMATIC) ? 128 : 0);


            DoubleBondStereochemistry stereo = stereochemistry[idx];

            if(stereo != null)
            {
                int idx0 = molecule.indexOf(bond.getAtom(0));
                int idx1 = molecule.indexOf(bond.getAtom(1));

                if(centers.elementType(idx0) != Type.Tricoordinate || centers.elementType(idx1) != Type.Tricoordinate)
                    continue;

                int flag = getDoubleBondStereoValue(bond, stereo.getBonds(), stereochemistry[idx].getStereo());
                array[offset + 3] |= flag;
            }
            else if(isDoubleBondStereochemistry(idx))
            {
                array[offset + 3] |= 0x03;
            }


            offset += bondBlockSize;
        }

        return array;
    }


    public byte[] hBondsAsBytes()
    {
        int count = 0;

        for(IBond bond : molecule.bonds())
        {
            if(bond.getAtom(0).getSymbol().equals("H") || bond.getAtom(1).getSymbol().equals("H"))
            {
                assert bond.getOrder() == Order.SINGLE || bond.getOrder() == Order.UNSET;
                assert !bond.getFlag(CDKConstants.ISAROMATIC);

                count++;
            }
        }


        byte[] array = new byte[count * hBondBlockSize];
        int offset = 0;

        for(IBond bond : molecule.bonds())
        {
            assert bond.getAtomCount() == 2;

            IAtom a1 = bond.getAtom(0);
            IAtom a2 = bond.getAtom(1);

            if(!a1.getSymbol().equals("H") && !a2.getSymbol().equals("H"))
                continue;

            int a1idx = molecule.indexOf(a1);
            int a2idx = molecule.indexOf(a2);
            assert a1idx != a2idx;

            if(bond.getOrder() == Order.SINGLE)
            {
                if(a1idx > a2idx)
                {
                    int tmp = a1idx;
                    a1idx = a2idx;
                    a2idx = tmp;
                }
            }
            else
            {
                if(a1idx < a2idx)
                {
                    int tmp = a1idx;
                    a1idx = a2idx;
                    a2idx = tmp;
                }
            }

            array[offset + 0] = (byte) (a1idx % 256);
            array[offset + 1] = (byte) (a1idx / 256 << 4 | a2idx / 256);
            array[offset + 2] = (byte) (a2idx % 256);

            offset += hBondBlockSize;
        }

        return array;
    }


    private boolean isTetrahedralChirality(int index)
    {
        if(!centers.isStereocenter(index))
            return false;

        if(centers.elementType(index) != Type.Tetracoordinate)
            return false;

        return true;
    }


    private boolean isExtendedTetrahedral(int index)
    {
        if(!centers.isStereocenter(index))
            return false;

        if(centers.elementType(index) != Type.Bicoordinate)
            return false;

        IAtom focus = molecule.getAtom(index);
        List<IBond> bonds = molecule.getConnectedBondsList(focus);

        if(bonds.size() != 2)
            return false;

        for(int i = 0; i < bonds.size(); i++)
            if(bonds.get(i).getOrder() != Order.DOUBLE)
                return false;


        /* find "left" terninal */
        IAtom leftPrevious = focus;
        IAtom leftTerminal = bonds.get(0).getOther(focus);
        int leftLength = 0;

        if(!centers.isStereocenter(molecule.indexOf(leftTerminal)))
            return false;

        while(true)
        {
            if(centers.elementType(molecule.indexOf(leftTerminal)) != Type.Bicoordinate)
                break;

            List<IAtom> neighbors = molecule.getConnectedAtomsList(leftTerminal);
            neighbors.remove(leftPrevious);

            if(neighbors.size() != 1)
                break;

            IAtom candidate = neighbors.get(0);

            if(molecule.getBond(leftTerminal, candidate).getOrder() != Order.DOUBLE)
                break;

            if(!centers.isStereocenter(molecule.indexOf(candidate)))
                break;

            leftPrevious = leftTerminal;
            leftTerminal = neighbors.get(0);
            leftLength++;
        }

        if(centers.elementType(molecule.indexOf(leftTerminal)) != Type.Tricoordinate)
            return false;


        /* find "right" terninal */
        IAtom rightPrevious = focus;
        IAtom rightTerminal = bonds.get(1).getOther(focus);
        int rightLength = 0;

        if(!centers.isStereocenter(molecule.indexOf(rightTerminal)))
            return false;

        while(true)
        {
            if(centers.elementType(molecule.indexOf(rightTerminal)) != Type.Bicoordinate)
                break;

            List<IAtom> neighbors = molecule.getConnectedAtomsList(rightTerminal);
            neighbors.remove(rightPrevious);

            if(neighbors.size() != 1)
                break;

            IAtom candidate = neighbors.get(0);

            if(molecule.getBond(rightTerminal, candidate).getOrder() != Order.DOUBLE)
                break;

            if(!centers.isStereocenter(molecule.indexOf(candidate)))
                break;

            rightPrevious = rightTerminal;
            rightTerminal = neighbors.get(0);
            rightLength++;
        }

        if(centers.elementType(molecule.indexOf(rightTerminal)) != Type.Tricoordinate)
            return false;


        return leftLength == rightLength;
    }


    private boolean isDoubleBondStereochemistry(int index)
    {
        IBond bond = molecule.getBond(index);

        if(bond.getOrder() != Order.DOUBLE)
            return false;

        for(IAtom atom : bond.atoms())
        {
            int idx = molecule.indexOf(atom);

            if(!centers.isStereocenter(idx))
                return false;

            if(centers.elementType(idx) != Type.Tricoordinate)
                return false;
        }

        return true;
    }


    private int getChiralityValue(IAtom center, IAtom[] ligands, Stereo stereo)
    {
        assert ligands.length == 4;

        ligands = ligands.clone();

        for(int i = 0; i < 4; i++)
            if(ligands[i].getSymbol().equals("H"))
                ligands[i] = center;


        int[] indexes = new int[4];

        for(int i = 0; i < 4; i++)
            indexes[i] = molecule.indexOf(ligands[i]);

        int order = 0;

        for(int i = 0; i < 4; i++)
        {
            int value = 0;

            for(int j = 0; j < 4; j++)
                if(indexes[j] <= indexes[i])
                    value++;

            order = (order << 4) + value;
        }

        boolean reverse = true;

        for(int i = 0; i < validReorder.length; i++)
            if(validReorder[i] == order)
                reverse = false;

        if(reverse)
            return stereo == Stereo.CLOCKWISE ? 2 : 1;
        else
            return stereo == Stereo.CLOCKWISE ? 1 : 2;
    }


    private int getDoubleBondStereoValue(IBond bond, IBond[] bonds, Conformation conformation)
    {
        if(!bonds[0].contains(bond.getAtom(0)))
        {
            bonds = bonds.clone();

            IBond tmp = bonds[0];
            bonds[0] = bonds[1];
            bonds[1] = tmp;
        }

        int[] ligands = new int[4];

        IAtom a0 = bonds[0].getOther(bond.getAtom(0));
        IAtom a1 = bonds[1].getOther(bond.getAtom(1));

        if(a0.getSymbol().equals("H"))
            a0 = bond.getAtom(0);

        if(a1.getSymbol().equals("H"))
            a1 = bond.getAtom(1);

        ligands[0] = molecule.indexOf(a0);
        ligands[1] = molecule.indexOf(bond.getAtom(0));
        ligands[2] = molecule.indexOf(a1);
        ligands[3] = molecule.indexOf(bond.getAtom(1));


        for(IBond b : molecule.bonds())
        {
            if(b == bond || b == bonds[0] || b == bonds[1])
                continue;

            if(b.contains(bond.getAtom(0)))
            {
                IAtom other = b.getOther(bond.getAtom(0));

                if(!other.getSymbol().equals("H"))
                    ligands[1] = molecule.indexOf(other);
            }

            if(b.contains(bond.getAtom(1)))
            {
                IAtom other = b.getOther(bond.getAtom(1));

                if(!other.getSymbol().equals("H"))
                    ligands[3] = molecule.indexOf(other);
            }
        }


        boolean reverse = false;

        if(ligands[0] > ligands[1])
            reverse = !reverse;

        if(ligands[2] > ligands[3])
            reverse = !reverse;

        if(reverse)
            return conformation == Conformation.OPPOSITE ? 2 : 1;
        else
            return conformation == Conformation.TOGETHER ? 1 : 2;
    }
}
