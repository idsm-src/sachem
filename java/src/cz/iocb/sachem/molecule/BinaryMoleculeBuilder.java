/*
 * Copyright (C) 2015-2020 Jakub Galgonek   galgonek@uochb.cas.cz
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
package cz.iocb.sachem.molecule;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.util.LinkedList;
import java.util.List;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry.Conformation;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;
import org.openscience.cdk.io.MDLV2000Writer.SPIN_MULTIPLICITY;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.sgroup.Sgroup;
import org.openscience.cdk.sgroup.SgroupKey;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.ExtendedCisTrans;
import org.openscience.cdk.stereo.ExtendedTetrahedral;
import org.openscience.cdk.stereo.TetrahedralChirality;
import cz.iocb.sachem.molecule.BinaryMolecule.SpecialRecordType;
import cz.iocb.sachem.molecule.BinaryMolecule.VariableLengthRecordType;
import cz.iocb.sachem.molecule.Molecule.AtomType;
import cz.iocb.sachem.molecule.Molecule.BondStereo;
import cz.iocb.sachem.molecule.Molecule.BondType;
import cz.iocb.sachem.molecule.Molecule.SgroupConnectivity;
import cz.iocb.sachem.molecule.Molecule.SgroupSubtype;
import cz.iocb.sachem.molecule.Molecule.SgroupType;
import cz.iocb.sachem.molecule.Molecule.TetrahedralStereo;



public class BinaryMoleculeBuilder
{
    public static String STEREO_PROPERTY = "STEREO";
    public static String IGNORE_STEREO = "IGNORE";
    public static String HYDROGEN_OFFSET = "HYDROGEN_OFFSET";
    private static int maxID = 4095;
    private static String BOND_NUMBER = "BOND_NUMBER";
    private static int[] validReorder = { 0x1234, 0x1423, 0x1342, 0x2314, 0x2431, 0x2143, 0x3124, 0x3412, 0x3241,
            0x4213, 0x4321, 0x4132 };


    public static byte[] asBytes(IAtomContainer molecule, boolean writeImplicitH) throws CDKException
    {
        @SuppressWarnings("rawtypes")
        IStereoElement[] tetrahedralChirality = new IStereoElement[molecule.getAtomCount()];

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


        @SuppressWarnings("rawtypes")
        IStereoElement[] doubleBondStereo = new IStereoElement[molecule.getBondCount()];

        for(@SuppressWarnings("rawtypes")
        IStereoElement element : molecule.stereoElements())
        {
            if(element instanceof DoubleBondStereochemistry)
            {
                DoubleBondStereochemistry e = (DoubleBondStereochemistry) element;
                doubleBondStereo[molecule.indexOf(e.getStereoBond())] = element;
            }
            else if(element instanceof ExtendedCisTrans)
            {
                ExtendedCisTrans e = (ExtendedCisTrans) element;
                doubleBondStereo[molecule.indexOf(e.getFocus())] = element;
            }
        }


        int xAtomCount = 0;
        int cAtomCount = 0;
        int hAtomCount = 0;
        int xBondCount = 0;
        int specialCount = 0;
        int sgroupCount = 0;

        for(IAtom a : molecule.atoms())
        {
            if(a instanceof IPseudoAtom)
                xAtomCount++;
            else if(a.getAtomicNumber() == Molecule.AtomType.H)
                hAtomCount++;
            else if(a.getAtomicNumber() == Molecule.AtomType.C)
                cAtomCount++;
            else
                xAtomCount++;

            if(writeImplicitH && a.getImplicitHydrogenCount() != null)
                hAtomCount += a.getImplicitHydrogenCount();

            if(a.getFormalCharge() != 0)
                specialCount++;

            if(a.getMassNumber() != null && a.getMassNumber() != -1)
                specialCount++;

            if(a.getProperty(CDKConstants.SPIN_MULTIPLICITY) != null)
                specialCount++;
        }

        for(IBond b : molecule.bonds())
        {
            boolean hBond = false;

            for(IAtom a : b.atoms())
                if(a.getAtomicNumber() == Molecule.AtomType.H && molecule.getConnectedBondsCount(a) == 1)
                    hBond = true;

            if(!hBond)
                xBondCount++;
        }

        for(IAtom atom : molecule.atoms())
            if(Boolean.TRUE.equals(atom.getProperty(MoleculeCreator.STEREO_FLAG)))
                specialCount++;

        for(IBond bond : molecule.bonds())
            if(Boolean.TRUE.equals(bond.getProperty(MoleculeCreator.STEREO_FLAG)))
                specialCount++;

        List<Sgroup> sgroups = molecule.getProperty(CDKConstants.CTAB_SGROUPS);

        if(sgroups != null)
            for(Sgroup sgroup : sgroups)
                if(getSgroupType(sgroup) != SgroupType.NONE)
                    sgroupCount++;

        if(sgroupCount > 0)
            specialCount++;


        if(xAtomCount + cAtomCount + hAtomCount > maxID || xBondCount + hAtomCount > maxID || specialCount > maxID)
            throw new CDKException("molecule is too big");


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
        List<IAtom> unknownAtoms = new LinkedList<IAtom>();

        for(int i = 0; i < xAtomCount; i++)
        {
            IAtom atom = molecule.getAtom(i);

            if(atom instanceof IPseudoAtom)
            {
                String label = ((IPseudoAtom) atom).getLabel();
                int type;

                if(label.length() == 1 && label.charAt(0) >= 'A' && label.charAt(0) <= 'Z')
                    type = -label.charAt(0);
                else if(label.matches("R[1-9#]?") || label.equals("*"))
                    type = AtomType.R;
                else if(label.equals("G*") || label.equals("G\\"))
                    type = AtomType.G;
                else if(label.equals("Ps"))
                    type = AtomType.POSITRONIUM;
                else if(label.equals("e"))
                    type = AtomType.ELECTRON;
                else if(label.equals("hv"))
                    type = AtomType.PHOTON;
                else if(label.equals("Mu"))
                    type = AtomType.MUONIUM;
                else if(label.equals("ACP"))
                    type = AtomType.ACP;
                else if(label.equals("Enz"))
                    type = AtomType.ENZYME;
                else if(label.equals(""))
                    type = AtomType.EMPTY;
                else
                    type = AtomType.UNKNOWN;

                stream.write(type);

                if(type == AtomType.UNKNOWN)
                    unknownAtoms.add(atom);
            }
            else
            {
                int num = molecule.getAtom(i).getAtomicNumber();
                assert num != Molecule.AtomType.C && num != Molecule.AtomType.H;

                assert num > 0 && num < 128;

                stream.write(num);
            }
        }


        int[] bondsCount = new int[molecule.getAtomCount()];

        /* write bonds */
        int heavyAtomCount = cAtomCount + xAtomCount;
        int bondNumber = 0;

        for(IBond b : molecule.bonds())
        {
            int a1idx = molecule.indexOf(b.getAtom(0));
            int a2idx = molecule.indexOf(b.getAtom(1));
            assert a1idx != a2idx;

            if(a1idx < heavyAtomCount && a2idx < heavyAtomCount)
            {
                b.setProperty(BOND_NUMBER, bondNumber++);
                bondsCount[a1idx]++;
                bondsCount[a2idx]++;

                stream.write(a1idx % 256);
                stream.write(a1idx / 256 << 4 | a2idx / 256);
                stream.write(a2idx % 256);
                stream.write(getBondType(b));
            }
        }


        /* write special bond number */
        for(IBond b : molecule.bonds())
        {
            IAtom a1 = b.getAtom(0);
            IAtom a2 = b.getAtom(1);

            int a1idx = molecule.indexOf(a1);
            int a2idx = molecule.indexOf(a2);
            assert a1idx != a2idx;

            if((a1idx >= heavyAtomCount || a2idx >= heavyAtomCount)
                    && (a1idx < heavyAtomCount || molecule.getConnectedBondsCount(a1) > 1)
                    && (a2idx < heavyAtomCount || molecule.getConnectedBondsCount(a2) > 1))
            {
                b.setProperty(BOND_NUMBER, bondNumber++);
                bondsCount[a1idx]++;
                bondsCount[a2idx]++;

                stream.write(a1idx % 256);
                stream.write(a1idx / 256 << 4 | a2idx / 256);
                stream.write(a2idx % 256);
                stream.write(getBondType(b));
            }
        }


        /* write hydrogen bonds */
        for(int i = xAtomCount + cAtomCount; i < molecule.getAtomCount(); i++)
        {
            IAtom atom = molecule.getAtom(i);

            List<IBond> list = molecule.getConnectedBondsList(atom);

            if(list.size() == 1)
            {
                IBond b = list.get(0);
                IAtom other = b.getOther(atom);
                int idx = molecule.indexOf(other);

                if(other.getAtomicNumber() == Molecule.AtomType.H && i > idx
                        && molecule.getConnectedBondsCount(other) == 1)
                {
                    // empty record
                    stream.write(0);
                    stream.write(0);
                }
                else
                {
                    b.setProperty(BOND_NUMBER, bondNumber++);
                    bondsCount[i]++;
                    bondsCount[idx]++;

                    stream.write(getBondType(molecule.getBond(atom, other)) << 4 | idx / 256);
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
        int offset = molecule.getAtomCount();

        if(writeImplicitH)
        {
            for(int idx = 0; idx < molecule.getAtomCount(); idx++)
            {
                IAtom atom = molecule.getAtom(idx);

                if(atom.getImplicitHydrogenCount() == null)
                    continue;

                atom.setProperty(HYDROGEN_OFFSET, offset);

                for(int h = 0; h < atom.getImplicitHydrogenCount(); h++)
                {
                    offset++;
                    bondsCount[idx]++;

                    stream.write(BondType.SINGLE << 4 | idx / 256);
                    stream.write(idx % 256);
                }
            }
        }


        for(int count : bondsCount)
            if(count > 16)
                throw new CDKException("too high atom valence");


        /* write charges, isotopes & radicals */
        for(IAtom a : molecule.atoms())
        {
            if(a.getFormalCharge() != 0)
            {
                int idx = molecule.indexOf(a);

                stream.write(SpecialRecordType.CHARGE << 4 | idx / 256);
                stream.write(idx % 256);
                stream.write(a.getFormalCharge());
            }

            if(a.getMassNumber() != null && a.getMassNumber() != -1)
            {
                int mass = a.getMassNumber() - a.getAtomicNumber() + 1;

                if(mass < 1)
                    throw new CDKException("wrong isotope number");

                if(mass > 256)
                    throw new CDKException("too high isotope number");

                int idx = molecule.indexOf(a);

                stream.write(SpecialRecordType.ISOTOPE << 4 | idx / 256);
                stream.write(idx % 256);
                stream.write(mass);
            }

            if(a.getProperty(CDKConstants.SPIN_MULTIPLICITY) != null)
            {
                SPIN_MULTIPLICITY multiplicity = (SPIN_MULTIPLICITY) a.getProperty(CDKConstants.SPIN_MULTIPLICITY);

                int idx = molecule.indexOf(a);

                stream.write(SpecialRecordType.RADICAL << 4 | idx / 256);
                stream.write(idx % 256);
                stream.write(multiplicity.getValue());
            }
        }


        /* write stereo centers */
        for(int idx = 0; idx < molecule.getAtomCount(); idx++)
        {
            if(Boolean.TRUE.equals(molecule.getAtom(idx).getProperty(MoleculeCreator.STEREO_FLAG)))
            {
                byte flag = TetrahedralStereo.UNDEFINED;

                if(tetrahedralChirality[idx] != null)
                {
                    if(tetrahedralChirality[idx] instanceof TetrahedralChirality)
                        flag = getTetrahedralChiralityType(molecule, (TetrahedralChirality) tetrahedralChirality[idx]);
                    else
                        flag = getExtendedTetrahedralType(molecule, (ExtendedTetrahedral) tetrahedralChirality[idx]);
                }

                stream.write(SpecialRecordType.TETRAHEDRAL_STEREO << 4 | idx / 256);
                stream.write(idx % 256);
                stream.write(flag);
            }
        }


        /* write double bond stereo */
        for(int idx = 0; idx < molecule.getBondCount(); idx++)
        {
            IBond bond = molecule.getBond(idx);

            if(Boolean.TRUE.equals(bond.getProperty(MoleculeCreator.STEREO_FLAG)))
            {
                byte flag = BondStereo.UNDEFINED;

                if(doubleBondStereo[idx] != null)
                {
                    if(doubleBondStereo[idx] instanceof DoubleBondStereochemistry)
                        flag = getDoubleBondStereochemistryType(molecule,
                                (DoubleBondStereochemistry) doubleBondStereo[idx]);
                    else
                        flag = getExtendedCisTransType(molecule, (ExtendedCisTrans) doubleBondStereo[idx]);
                }

                int index = bond.getProperty(BOND_NUMBER);
                stream.write(SpecialRecordType.BOND_STEREO << 4 | index / 256);
                stream.write(index % 256);
                stream.write(flag);
            }
        }


        /* write SGroup count */
        if(sgroups != null)
        {
            if(sgroupCount > maxID)
                throw new CDKException("molecule is too big");

            stream.write(SpecialRecordType.SGROUP_COUNT << 4 | sgroupCount / 256);
            stream.write(sgroupCount % 256);
            stream.write(0);
        }


        /* write pseudo atom labels */
        for(IAtom atom : unknownAtoms)
        {
            String label = ((IPseudoAtom) atom).getLabel();

            int idx = molecule.indexOf(atom);
            byte[] bytes = label.getBytes(StandardCharsets.UTF_8);
            int size = 7 + bytes.length;

            stream.write(size / (256 * 256 * 256));
            stream.write(size / (256 * 256) % 256);
            stream.write(size / 256 % 256);
            stream.write(size % 256);

            stream.write(VariableLengthRecordType.ATOM_LABEL);
            stream.write(idx / 256);
            stream.write(idx % 256);

            for(byte b : bytes)
                stream.write(b);
        }


        /* write SGroups */
        if(sgroups != null)
        {
            for(Sgroup sgroup : sgroups)
            {
                if(getSgroupType(sgroup) == SgroupType.NONE)
                    continue;

                int atomCount = sgroup.getAtoms().size();
                int bondCount = sgroup.getBonds().size();

                if(writeImplicitH)
                    for(IAtom atom : sgroup.getAtoms())
                        if(atom.getImplicitHydrogenCount() != null)
                            atomCount += atom.getImplicitHydrogenCount();

                int size = 12 + 2 * atomCount + 4 * bondCount;

                stream.write(size / (256 * 256 * 256));
                stream.write(size / (256 * 256) % 256);
                stream.write(size / 256 % 256);
                stream.write(size % 256);

                stream.write(VariableLengthRecordType.SGROUP);
                stream.write(getSgroupType(sgroup));
                stream.write(getSgroupSubtype(sgroup));
                stream.write(getSgroupConnectivity(sgroup));

                stream.write(atomCount / 256);
                stream.write(atomCount % 256);
                stream.write(bondCount / 256);
                stream.write(bondCount % 256);

                for(IAtom atom : sgroup.getAtoms())
                {
                    int index = molecule.indexOf(atom);

                    stream.write(index / 256);
                    stream.write(index % 256);

                    if(writeImplicitH)
                    {
                        if(atom.getImplicitHydrogenCount() != null)
                        {
                            int base = atom.getProperty(HYDROGEN_OFFSET);

                            for(int h = 0; h < atom.getImplicitHydrogenCount(); h++)
                            {
                                stream.write((base + h) / 256);
                                stream.write((base + h) % 256);
                            }
                        }
                    }
                }

                for(IBond bond : sgroup.getBonds())
                {
                    for(IAtom atom : bond.atoms())
                    {
                        int index = molecule.indexOf(atom);
                        stream.write(index / 256);
                        stream.write(index % 256);
                    }
                }
            }
        }


        return stream.toByteArray();
    }


    public static byte getBondType(IBond bond) throws CDKException
    {
        if(bond.isAromatic())
            return BondType.AROMATIC;

        if(bond instanceof QueryBond)
        {
            switch(((QueryBond) bond).getExpression().type())
            {
                case SINGLE_OR_DOUBLE:
                    return BondType.SINGLE_OR_DOUBLE;
                case SINGLE_OR_AROMATIC:
                    return BondType.SINGLE_OR_AROMATIC;
                case DOUBLE_OR_AROMATIC:
                    return BondType.DOUBLE_OR_AROMATIC;
                default:
                    return BondType.ANY;
            }
        }

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
                return BondType.ANY;
        }

        throw new CDKException("unknown bond type");
    }


    private static byte getTetrahedralChiralityType(IAtomContainer molecule, TetrahedralChirality chirality)
    {
        IAtom center = chirality.getChiralAtom();
        IAtom[] ligands = chirality.getLigands();
        Stereo stereo = chirality.getStereo();

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
            return stereo == Stereo.CLOCKWISE ? TetrahedralStereo.ANTI_CLOCKWISE : TetrahedralStereo.CLOCKWISE;
        else
            return stereo == Stereo.CLOCKWISE ? TetrahedralStereo.CLOCKWISE : TetrahedralStereo.ANTI_CLOCKWISE;
    }


    private static byte getExtendedTetrahedralType(IAtomContainer molecule, ExtendedTetrahedral chirality)
    {
        IAtom[] terminals = chirality.findTerminalAtoms(molecule);
        IAtom[] ligands = chirality.peripherals();
        Stereo stereo = chirality.winding();

        assert ligands.length == 4;
        int[] indexes = new int[4];

        for(int i = 0; i < 4; i++)
        {
            if(ligands[i] != terminals[0] && ligands[i] != terminals[1])
                indexes[i] = molecule.indexOf(ligands[i]);
            else
                indexes[i] = Integer.MAX_VALUE;
        }

        boolean reverse = true;

        if(indexes[0] > indexes[1])
            reverse = !reverse;

        if(indexes[2] > indexes[3])
            reverse = !reverse;

        if(reverse)
            return stereo == Stereo.CLOCKWISE ? TetrahedralStereo.ANTI_CLOCKWISE : TetrahedralStereo.CLOCKWISE;
        else
            return stereo == Stereo.CLOCKWISE ? TetrahedralStereo.CLOCKWISE : TetrahedralStereo.ANTI_CLOCKWISE;
    }


    private static byte getDoubleBondStereochemistryType(IAtomContainer molecule, DoubleBondStereochemistry stereo)
    {
        IBond bond = stereo.getStereoBond();
        IBond[] bonds = stereo.getBonds();
        Conformation conformation = stereo.getStereo();

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
            return conformation == Conformation.TOGETHER ? BondStereo.OPPOSITE : BondStereo.TOGETHER;
        else
            return conformation == Conformation.TOGETHER ? BondStereo.TOGETHER : BondStereo.OPPOSITE;
    }


    private static byte getExtendedCisTransType(IAtomContainer molecule, ExtendedCisTrans stereo)
    {
        List<IBond> carriers = stereo.getCarriers();
        int conformation = stereo.getConfigOrder();


        int[] ligands = new int[4];

        for(int i = 0; i < 2; i++)
        {
            IBond bond = stereo.getFocus();
            IAtom terminal = bond.getAtom(i);

            while(true)
            {
                List<IBond> bonds = molecule.getConnectedBondsList(terminal);

                if(bonds.size() != 2)
                    break;

                IBond next = bonds.get(0) == bond ? bonds.get(1) : bonds.get(0);

                if(next.getOrder() != IBond.Order.DOUBLE)
                    break;

                terminal = next.getOther(terminal);
                bond = next;
            }

            IBond carrier = carriers.get(0).contains(terminal) ? carriers.get(0) : carriers.get(1);
            IBond cocarier = null;

            for(IBond b : molecule.getConnectedBondsList(terminal))
                if(b != carrier && b != bond)
                    cocarier = b;

            ligands[2 * i] = molecule.indexOf(carrier.getOther(terminal));
            ligands[2 * i + 1] = cocarier != null ? molecule.indexOf(cocarier.getOther(terminal)) : Integer.MAX_VALUE;
        }


        boolean reverse = false;

        if(ligands[0] > ligands[1])
            reverse = !reverse;

        if(ligands[2] > ligands[3])
            reverse = !reverse;


        if(reverse)
            return conformation == IStereoElement.TOGETHER ? BondStereo.OPPOSITE : BondStereo.TOGETHER;
        else
            return conformation == IStereoElement.TOGETHER ? BondStereo.TOGETHER : BondStereo.OPPOSITE;
    }


    private static byte getSgroupType(Sgroup sgroup)
    {
        switch(sgroup.getType())
        {
            case CtabStructureRepeatUnit:
                return SgroupType.SRU;

            case CtabModified:
                return SgroupType.MOD;

            case CtabMonomer:
                return SgroupType.MON;

            case CtabCopolymer:
                return SgroupType.COP;

            case CtabGeneric:
                return SgroupType.GEN;


            case CtabAnyPolymer:
                return SgroupType.ANY;

            case CtabCrossLink:
                return SgroupType.CRO;

            case CtabMer:
                return SgroupType.MER;

            case CtabGraft:
                return SgroupType.GRA;


            case CtabComponent:
                return SgroupType.COM;

            case CtabFormulation:
                return SgroupType.FOR;

            case CtabMixture:
                return SgroupType.MIX;

            default:
                return SgroupType.NONE;
        }
    }


    private static byte getSgroupSubtype(Sgroup sgroup)
    {
        if(sgroup == null)
            return SgroupSubtype.NONE;

        String subtype = sgroup.getValue(SgroupKey.CtabSubType);

        if(subtype == null)
            return SgroupSubtype.NONE;

        switch(subtype)
        {
            case "ALT":
                return SgroupSubtype.ALT;

            case "RAN":
                return SgroupSubtype.RAN;

            case "BLO":
                return SgroupSubtype.BLO;

            default:
                return SgroupSubtype.UNKNOWN;
        }
    }


    private static byte getSgroupConnectivity(Sgroup sgroup)
    {
        if(sgroup == null)
            return SgroupConnectivity.NONE;

        String connectivity = sgroup.getValue(SgroupKey.CtabConnectivity);

        if(connectivity == null)
            return SgroupConnectivity.NONE;

        switch(connectivity)
        {
            case "HH":
                return SgroupConnectivity.HH;

            case "HT":
                return SgroupConnectivity.HT;

            case "EU":
                return SgroupConnectivity.EU;

            default:
                return SgroupConnectivity.UNKNOWN;
        }
    }
}
