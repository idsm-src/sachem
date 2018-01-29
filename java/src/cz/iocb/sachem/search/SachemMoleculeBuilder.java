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
package cz.iocb.sachem.search;

import java.io.ByteArrayOutputStream;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry.Conformation;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;
import org.openscience.cdk.isomorphism.matchers.CTFileQueryBond;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.ExtendedTetrahedral;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.stereo.Stereocenters.Type;
import org.openscience.cdk.stereo.TetrahedralChirality;
import org.openscience.cdk.tools.periodictable.PeriodicTable;



public class SachemMoleculeBuilder
{
    public static enum BondType
    {
        NONE(0),
        SINGLE(1),
        DOUBLE(2),
        TRIPLE(3),
        QUADRUPLE(4),
        QUINTUPLE(5),
        SEXTUPLE(6),
        AROMATIC(11),
        SINGLE_OR_DOUBLE(12),
        SINGLE_OR_AROMATIC(13),
        DOUBLE_OR_AROMATIC(14),
        ANY(15);

        private int value;

        private BondType(int value)
        {
            this.value = value;
        };

        public int getValue()
        {
            return value;
        }
    }


    static enum SpecialRecordType
    {
        CHARGE(0), ISOTOPE(1), TETRAHEDRAL_STEREO(2), BOND_STEREO(3),;

        private int value;

        private SpecialRecordType(int value)
        {
            this.value = value;
        };

        public int getValue()
        {
            return value;
        }
    }


    static enum TetrahedralStereoType
    {
        NONE(0), CLOCKWISE(1), ANTI_CLOCKWISE(2), UNDEFINED(3);

        private int value;

        private TetrahedralStereoType(int value)
        {
            this.value = value;
        };

        public int getValue()
        {
            return value;
        }
    }


    static enum BondStereoType
    {
        NONE(0), OPPOSITE(1), TOGETHER(2), UNDEFINED(3);

        private int value;

        private BondStereoType(int value)
        {
            this.value = value;
        };

        public int getValue()
        {
            return value;
        }
    }


    private static String BOND_NUMBER = "BOND_NUMBER";
    private static int[] validReorder = { 0x1234, 0x1423, 0x1342, 0x2314, 0x2431, 0x2143, 0x3124, 0x3412, 0x3241,
            0x4213, 0x4321, 0x4132 };

    private IAtomContainer molecule;
    private Stereocenters centers;
    @SuppressWarnings("rawtypes") private IStereoElement[] tetrahedralChirality;
    private DoubleBondStereochemistry[] doubleBondStereo;


    public SachemMoleculeBuilder(IAtomContainer molecule) throws CDKException
    {
        this.molecule = molecule;

        try
        {
            centers = Stereocenters.of(molecule);
        }
        catch(Throwable e)
        {
        }


        tetrahedralChirality = new IStereoElement[molecule.getAtomCount()];

        for(@SuppressWarnings("rawtypes")
        IStereoElement element : molecule.stereoElements())
        {
            if(element instanceof TetrahedralChirality)
            {
                int index = molecule.indexOf(((TetrahedralChirality) element).getChiralAtom());
                tetrahedralChirality[index] = element;
            }
            else if(element instanceof ExtendedTetrahedral)
            {
                int index = molecule.indexOf(((ExtendedTetrahedral) element).focus());
                tetrahedralChirality[index] = element;
            }
        }


        doubleBondStereo = new DoubleBondStereochemistry[molecule.getBondCount()];

        for(@SuppressWarnings("rawtypes")
        IStereoElement element : molecule.stereoElements())
        {
            if(element instanceof DoubleBondStereochemistry)
            {
                DoubleBondStereochemistry e = (DoubleBondStereochemistry) element;
                doubleBondStereo[molecule.indexOf(e.getStereoBond())] = e;
            }
        }
    }


    public byte[] asBytes(boolean writeImplicitH) throws CDKException
    {
        int xAtomCount = 0;
        int cAtomCount = 0;
        int hAtomCount = 0;
        int xBondCount = 0;
        int specialCount = 0;

        for(IAtom a : molecule.atoms())
        {
            if(a instanceof IPseudoAtom)
                xAtomCount++;
            else if(a.getSymbol().equals("H"))
                hAtomCount++;
            else if(a.getSymbol().equals("C"))
                cAtomCount++;
            else
                xAtomCount++;

            if(writeImplicitH && a.getImplicitHydrogenCount() != null)
                hAtomCount += a.getImplicitHydrogenCount();

            if(a.getFormalCharge() != 0)
                specialCount++;

            if(a.getMassNumber() != null)
                specialCount++;
        }

        for(IBond b : molecule.bonds())
        {
            boolean hBond = false;

            for(IAtom a : b.atoms())
                if(a.getSymbol().equals("H") && molecule.getConnectedBondsCount(a) == 1)
                    hBond = true;

            if(!hBond)
                xBondCount++;
        }

        for(int idx = 0; idx < molecule.getAtomCount(); idx++)
            if(isTetrahedralChirality(idx) || isExtendedTetrahedral(idx))
                specialCount++;

        for(int idx = 0; idx < molecule.getBondCount(); idx++)
            if(isDoubleBondStereochemistry(idx))
                specialCount++;


        int length = 5 * 2 + xAtomCount + 4 * xBondCount + 2 * hAtomCount + 3 * specialCount;
        ByteArrayOutputStream stream = new ByteArrayOutputStream(length);

        stream.write(xAtomCount >> 8);
        stream.write(xAtomCount & 0xFF);

        stream.write(cAtomCount >> 8);
        stream.write(cAtomCount & 0xFF);

        stream.write(hAtomCount >> 8);
        stream.write(hAtomCount & 0xFF);

        stream.write(xBondCount >> 8);
        stream.write(xBondCount & 0xFF);

        stream.write(specialCount >> 8);
        stream.write(specialCount & 0xFF);


        /* write atoms */
        for(int i = 0; i < xAtomCount; i++)
        {
            IAtom atom = molecule.getAtom(i);

            if(atom instanceof IPseudoAtom)
            {
                String label = ((IPseudoAtom) atom).getLabel();
                assert label.length() == 0 || label.length() > 1 || label.charAt(0) >= 0 && label.charAt(0) < 128;

                if(label.length() == 1)
                    stream.write(-label.charAt(0));
                else if(label.matches("R[1-9#]?"))
                    stream.write(-'R');
                else
                    stream.write(-'?');
            }
            else
            {
                String symbol = molecule.getAtom(i).getSymbol();
                assert !symbol.equals("C") && !symbol.equals("H");

                int num = PeriodicTable.getAtomicNumber(symbol);
                assert num > 0 && num < 128;

                stream.write(num);
            }
        }


        /* write bonds */
        int heavyAtomCount = cAtomCount + xAtomCount;
        int bondNumber = 0;

        for(IBond b : molecule.bonds())
        {
            IAtom a1 = b.getAtom(0);
            IAtom a2 = b.getAtom(1);

            int a1idx = molecule.indexOf(b.getAtom(0));
            int a2idx = molecule.indexOf(b.getAtom(1));
            assert a1idx != a2idx;

            if((a1idx < heavyAtomCount || molecule.getConnectedBondsCount(a1) > 1)
                    && (a2idx < heavyAtomCount || molecule.getConnectedBondsCount(a2) > 1))
            {
                b.setProperty(BOND_NUMBER, bondNumber++);

                stream.write(a1idx % 256);
                stream.write(a1idx / 256 << 4 | a2idx / 256);
                stream.write(a2idx % 256);
                stream.write(getBondType(b).getValue());
            }
        }


        /* write hydrogen bonds */
        for(int i = xAtomCount + cAtomCount; i < molecule.getAtomCount(); i++)
        {
            IAtom atom = molecule.getAtom(i);

            List<IBond> list = molecule.getConnectedBondsList(atom);
            assert list.size() <= 1;

            if(list.size() == 1)
            {
                IAtom other = list.get(0).getOther(atom);
                int idx = molecule.indexOf(other);

                if(other.getSymbol().equals("H") && i > idx)
                {
                    // empty record
                    stream.write(0);
                    stream.write(0);
                }
                else
                {
                    stream.write(getBondType(molecule.getBond(atom, other)).getValue() << 4 | idx / 256);
                    stream.write(idx % 256);
                }
            }
            else
            {
                // empty record
                stream.write(0);
                stream.write(0);
            }
        }


        /* write implicit hydrogen bonds */
        if(writeImplicitH)
        {
            for(int idx = 0; idx < molecule.getAtomCount(); idx++)
            {
                IAtom atom = molecule.getAtom(idx);

                if(atom.getImplicitHydrogenCount() == null)
                    continue;

                for(int h = 0; h < atom.getImplicitHydrogenCount(); h++)
                {
                    stream.write(BondType.SINGLE.getValue() << 4 | idx / 256);
                    stream.write(idx % 256);
                }
            }
        }


        /* write charges & isotopes */
        for(IAtom a : molecule.atoms())
        {
            if(a.getFormalCharge() != 0)
            {
                int idx = molecule.indexOf(a);

                stream.write(SpecialRecordType.CHARGE.getValue() << 4 | idx / 256);
                stream.write(idx % 256);
                stream.write(a.getFormalCharge());
            }

            if(a.getMassNumber() != null)
            {
                int mass = a.getMassNumber() - a.getAtomicNumber() + 1;

                if(mass < 1)
                    throw new CDKException("wrong isotope number");

                if(mass > 256)
                    throw new CDKException("too high isotope number");

                int idx = molecule.indexOf(a);

                stream.write(SpecialRecordType.ISOTOPE.getValue() << 4 | idx / 256);
                stream.write(idx % 256);
                stream.write(mass);
            }
        }


        /* write stereo centers */
        for(int idx = 0; idx < molecule.getAtomCount(); idx++)
        {
            if(isTetrahedralChirality(idx) || isExtendedTetrahedral(idx))
            {
                TetrahedralStereoType flag = TetrahedralStereoType.UNDEFINED;

                @SuppressWarnings("rawtypes")
                IStereoElement element = tetrahedralChirality[idx];

                if(element != null && element instanceof TetrahedralChirality)
                {
                    TetrahedralChirality chirality = (TetrahedralChirality) element;
                    flag = getChiralityValue(chirality.getChiralAtom(), chirality.getLigands(), chirality.getStereo());
                }
                else if(element != null && element instanceof ExtendedTetrahedral)
                {
                    ExtendedTetrahedral chirality = (ExtendedTetrahedral) element;
                    flag = getChiralityValue(chirality.focus(), chirality.peripherals(), chirality.winding());
                }

                stream.write(SpecialRecordType.TETRAHEDRAL_STEREO.getValue() << 4 | idx / 256);
                stream.write(idx % 256);
                stream.write(flag.getValue());
            }
        }


        /* write double bond stereo */
        for(int idx = 0; idx < molecule.getBondCount(); idx++)
        {
            IBond bond = molecule.getBond(idx);
            assert bond.getAtomCount() == 2;

            if(isDoubleBondStereochemistry(idx))
            {
                BondStereoType flag = BondStereoType.UNDEFINED;

                DoubleBondStereochemistry stereo = doubleBondStereo[idx];

                if(stereo != null)
                    flag = getDoubleBondStereoType(bond, stereo.getBonds(), doubleBondStereo[idx].getStereo());

                int index = bond.getProperty(BOND_NUMBER);
                stream.write(SpecialRecordType.BOND_STEREO.getValue() << 4 | index / 256);
                stream.write(index % 256);
                stream.write(flag.getValue());
            }
        }


        byte[] array = stream.toByteArray();
        assert array.length == length;

        return array;
    }


    public static BondType getBondType(IBond bond) throws CDKException
    {
        if(bond.isAromatic())
            return BondType.AROMATIC;

        switch(bond.getOrder())
        {
            case SINGLE:
                return BondType.SINGLE;
            case DOUBLE:
                return BondType.DOUBLE;
            case TRIPLE:
                return BondType.TRIPLE;
            case QUADRUPLE:
                return BondType.QUADRUPLE;
            case QUINTUPLE:
                return BondType.QUINTUPLE;
            case SEXTUPLE:
                return BondType.SEXTUPLE;
            case UNSET:
                if(bond instanceof CTFileQueryBond)
                {
                    switch(((CTFileQueryBond) bond).getType())
                    {
                        case SINGLE:
                            return BondType.SINGLE;
                        case DOUBLE:
                            return BondType.DOUBLE;
                        case TRIPLE:
                            return BondType.TRIPLE;
                        case AROMATIC:
                            return BondType.AROMATIC;
                        case SINGLE_OR_DOUBLE:
                            return BondType.SINGLE_OR_DOUBLE;
                        case SINGLE_OR_AROMATIC:
                            return BondType.SINGLE_OR_AROMATIC;
                        case DOUBLE_OR_AROMATIC:
                            return BondType.DOUBLE_OR_AROMATIC;
                        case ANY:
                            return BondType.ANY;
                    }
                }
                else
                {
                    return BondType.ANY;
                }
        }

        throw new CDKException("unknown bond type");
    }


    private boolean isTetrahedralChirality(int index)
    {
        if(tetrahedralChirality[index] != null)
            return true;

        if(centers == null)
            return false;

        if(!centers.isStereocenter(index))
            return false;

        if(centers.elementType(index) != Type.Tetracoordinate)
            return false;

        return true;
    }


    private boolean isExtendedTetrahedral(int index)
    {
        if(tetrahedralChirality[index] != null)
            return true;

        if(centers == null)
            return false;

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

        if(bond.isAromatic())
            return false;

        if(doubleBondStereo[index] != null)
            return true;

        if(centers == null)
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


    private TetrahedralStereoType getChiralityValue(IAtom center, IAtom[] ligands, Stereo stereo)
    {
        assert ligands.length == 4;
        int[] indexes = new int[4];

        for(int i = 0; i < 4; i++)
        {
            if(ligands[i] != center)
                indexes[i] = molecule.indexOf(ligands[i]);
            else
                indexes[i] = Integer.MAX_VALUE;
        }

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
            return stereo == Stereo.CLOCKWISE ? TetrahedralStereoType.ANTI_CLOCKWISE : TetrahedralStereoType.CLOCKWISE;
        else
            return stereo == Stereo.CLOCKWISE ? TetrahedralStereoType.CLOCKWISE : TetrahedralStereoType.ANTI_CLOCKWISE;
    }


    private BondStereoType getDoubleBondStereoType(IBond bond, IBond[] bonds, Conformation conformation)
    {
        if(!bonds[0].contains(bond.getAtom(0)))
            bonds = new IBond[] { bonds[1], bonds[0] };

        int[] ligands = new int[4];
        ligands[0] = molecule.indexOf(bonds[0].getOther(bond.getAtom(0)));
        ligands[1] = Integer.MAX_VALUE;
        ligands[2] = molecule.indexOf(bonds[1].getOther(bond.getAtom(1)));
        ligands[3] = Integer.MAX_VALUE;


        for(IBond b : molecule.bonds())
        {
            if(b == bond || b == bonds[0] || b == bonds[1])
                continue;

            if(b.contains(bond.getAtom(0)))
                ligands[1] = molecule.indexOf(b.getOther(bond.getAtom(0)));

            if(b.contains(bond.getAtom(1)))
                ligands[3] = molecule.indexOf(b.getOther(bond.getAtom(1)));
        }


        boolean reverse = false;

        if(ligands[0] > ligands[1])
            reverse = !reverse;

        if(ligands[2] > ligands[3])
            reverse = !reverse;

        if(reverse)
            return conformation == Conformation.TOGETHER ? BondStereoType.OPPOSITE : BondStereoType.TOGETHER;
        else
            return conformation == Conformation.TOGETHER ? BondStereoType.TOGETHER : BondStereoType.OPPOSITE;
    }
}
