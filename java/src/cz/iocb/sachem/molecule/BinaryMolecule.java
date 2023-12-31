package cz.iocb.sachem.molecule;

import java.util.ArrayList;



public final class BinaryMolecule extends Molecule
{
    public static abstract class SpecialRecordType
    {
        public static final byte CHARGE = 0;
        public static final byte ISOTOPE = 1;
        public static final byte TETRAHEDRAL_STEREO = 2;
        public static final byte BOND_STEREO = 3;
        public static final byte RADICAL = 4;
        public static final byte SGROUP_COUNT = 5;
    }


    public static abstract class VariableLengthRecordType
    {
        public static final byte SGROUP = 1;
        public static final byte ATOM_LABEL = 2;
    }


    private static final int BOND_BLOCK_SIZE = 4;
    private static final int HBOND_BLOCK_SIZE = 2;
    private static final int SPECIAL_BLOCK_SIZE = 3;

    private final int originalAtomCount;
    private final int originalBondCount;
    private final int atomCount;
    private final int bondCount;

    private final byte[] atomNumbers;
    private final byte[] bondTypes;
    private final int[][] bondLists;
    private final int[] contains;
    private final int[] bondMatrix;
    private final byte[] atomHydrogens;
    private final byte[] atomCharges;
    private final byte[] atomMasses;
    private final byte[] atomRadicalTypes;
    private final byte[] atomStereo;
    private final byte[] bondStereo;
    private final boolean[] restH;
    private final SGroup[] sgroups;
    private final AtomLabel[] labels;


    public BinaryMolecule(byte[] data, boolean[] restH, boolean extended, boolean withCharges, boolean withIsotopes,
            boolean withRadicals, boolean withStereo, boolean withSGroup, boolean ignoreChargedHydrogens,
            boolean ignoreHydrogenIsotopes, boolean ignoreHydrogenRadicals)
    {
        ignoreChargedHydrogens &= !extended;
        ignoreHydrogenIsotopes &= !extended;
        ignoreHydrogenRadicals &= !extended;

        int possition = 0;

        int xAtomCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);
        int cAtomCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);
        int hAtomCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);
        int xBondCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);
        int specialCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);

        int heavyAtomCount = xAtomCount + cAtomCount;
        int originalAtomCount = heavyAtomCount + hAtomCount;
        int originalBondCount = xBondCount + hAtomCount;
        int atomCount = heavyAtomCount;
        int bondCount = xBondCount;
        int sgroupCount = 0;
        int labelCount = 0;

        if(extended)
        {
            atomCount += hAtomCount;
            bondCount += hAtomCount;
        }

        byte[] atomNumbers = new byte[atomCount];
        byte[] atomHydrogens = new byte[atomCount];
        byte[] bondTypes = new byte[bondCount];
        byte[] atomCharges = withCharges ? new byte[atomCount] : null;
        byte[] atomMasses = withIsotopes ? new byte[atomCount] : null;
        byte[] atomRadicalTypes = withRadicals ? new byte[atomCount] : null;
        byte[] atomStereo = withStereo ? new byte[atomCount] : null;
        byte[] bondStereo = withStereo ? new byte[bondCount] : null;
        int[] contains = new int[bondCount * 2];
        int[] bondMatrix = new int[atomCount * atomCount];


        for(int i = 0; i < xAtomCount; i++)
            atomNumbers[i] = data[possition++];

        for(int i = xAtomCount; i < heavyAtomCount; i++)
            atomNumbers[i] = AtomType.C;

        for(int i = heavyAtomCount; i < atomCount; i++)
            atomNumbers[i] = AtomType.H;

        for(int i = 0; i < xAtomCount; i++)
            if(atomNumbers[i] == AtomType.UNKNOWN)
                labelCount++;


        for(int i = 0; i < atomCount * atomCount; i++)
            bondMatrix[i] = -1;

        int boundIdx = 0;

        ArrayList<ArrayList<Integer>> tmpBondLists = new ArrayList<ArrayList<Integer>>(atomCount);

        for(int i = 0; i < atomCount; i++)
            tmpBondLists.add(new ArrayList<Integer>(16));

        for(int i = 0; i < xBondCount; i++)
        {
            int b0 = Byte.toUnsignedInt(data[possition++]);
            int b1 = Byte.toUnsignedInt(data[possition++]);
            int b2 = Byte.toUnsignedInt(data[possition++]);
            bondTypes[i] = data[possition++];


            int x = b0 | (b1 << 4 & 0xF00);
            int y = b2 | (b1 << 8 & 0xF00);

            if(x >= heavyAtomCount && y < atomCount)
                atomHydrogens[y]++;

            if(y >= heavyAtomCount && x < atomCount)
                atomHydrogens[x]++;

            if(x >= atomCount || y >= atomCount)
            {
                bondCount--;
                continue;
            }

            tmpBondLists.get(x).add(y);
            tmpBondLists.get(y).add(x);

            bondMatrix[x * atomCount + y] = boundIdx;
            bondMatrix[y * atomCount + x] = boundIdx;

            contains[boundIdx * 2 + 0] = x;
            contains[boundIdx * 2 + 1] = y;

            boundIdx++;
        }


        boolean[] ignoredHydrogen = null;

        if(ignoreChargedHydrogens || ignoreHydrogenIsotopes || ignoreHydrogenRadicals)
        {
            int base = possition + hAtomCount * HBOND_BLOCK_SIZE;
            ignoredHydrogen = new boolean[hAtomCount];

            for(int i = 0; i < specialCount; i++)
            {
                int offset = base + i * SPECIAL_BLOCK_SIZE;
                int value = Byte.toUnsignedInt(data[offset + 0]) * 256 | Byte.toUnsignedInt(data[offset + 1]);
                int idx = value & 0xFFF;

                switch(Byte.toUnsignedInt(data[offset]) >> 4)
                {
                    case SpecialRecordType.CHARGE:
                        if(ignoreChargedHydrogens && idx >= heavyAtomCount)
                            ignoredHydrogen[idx - heavyAtomCount] = true;
                        break;

                    case SpecialRecordType.ISOTOPE:
                        if(ignoreHydrogenIsotopes && idx >= heavyAtomCount)
                            ignoredHydrogen[idx - heavyAtomCount] = true;
                        break;

                    case SpecialRecordType.RADICAL:
                        if(ignoreHydrogenRadicals && idx >= heavyAtomCount)
                            ignoredHydrogen[idx - heavyAtomCount] = true;
                        break;
                }
            }
        }


        for(int i = 0; i < hAtomCount; i++)
        {
            int offset = possition + i * HBOND_BLOCK_SIZE;
            int value = Byte.toUnsignedInt(data[offset + 0]) * 256 | Byte.toUnsignedInt(data[offset + 1]);

            if(value == 0)
            {
                originalBondCount--;

                if(extended)
                    bondCount--;

                continue;
            }


            int x = value & 0xFFF;

            if(x < atomCount && (!(ignoreChargedHydrogens || ignoreHydrogenIsotopes || ignoreHydrogenRadicals)
                    || ignoredHydrogen[i] == false))
                atomHydrogens[x]++;


            if(extended)
            {
                int y = heavyAtomCount + i;

                if(x >= heavyAtomCount)
                    atomHydrogens[y]++;

                bondTypes[boundIdx] = (byte) (Byte.toUnsignedInt(data[offset]) >> 4);


                tmpBondLists.get(x).add(y);
                tmpBondLists.get(y).add(x);

                bondMatrix[x * atomCount + y] = boundIdx;
                bondMatrix[y * atomCount + x] = boundIdx;

                contains[boundIdx * 2 + 0] = x;
                contains[boundIdx * 2 + 1] = y;

                boundIdx++;
            }
        }

        possition += hAtomCount * HBOND_BLOCK_SIZE;


        for(int i = 0; i < specialCount; i++)
        {
            int offset = possition + i * SPECIAL_BLOCK_SIZE;
            int value = Byte.toUnsignedInt(data[offset + 0]) * 256 | Byte.toUnsignedInt(data[offset + 1]);
            int idx = value & 0xFFF;

            switch(Byte.toUnsignedInt(data[offset]) >> 4)
            {
                case SpecialRecordType.CHARGE:
                    if(withCharges && idx < atomCount)
                        atomCharges[idx] = data[offset + 2];
                    break;

                case SpecialRecordType.ISOTOPE:
                    if(withIsotopes && idx < atomCount)
                        atomMasses[idx] = data[offset + 2];
                    break;

                case SpecialRecordType.RADICAL:
                    if(withRadicals && idx < atomCount)
                        atomRadicalTypes[idx] = data[offset + 2];
                    break;

                case SpecialRecordType.TETRAHEDRAL_STEREO:
                    if(withStereo && idx < atomCount)
                        atomStereo[idx] = data[offset + 2];
                    break;

                case SpecialRecordType.BOND_STEREO:
                    if(withStereo && idx < xBondCount) // not for H bond
                        bondStereo[idx] = data[offset + 2];
                    break;

                case SpecialRecordType.SGROUP_COUNT:
                    if(withSGroup)
                        sgroupCount = idx;
                    break;
            }
        }


        possition += specialCount * SPECIAL_BLOCK_SIZE;


        AtomLabel[] labels = null;

        if(labelCount > 0)
            labels = new AtomLabel[labelCount];

        for(int i = 0; i < labelCount; i++)
        {
            int size = Byte.toUnsignedInt(data[possition + 0]) << 24 | Byte.toUnsignedInt(data[possition + 1]) << 16
                    | Byte.toUnsignedInt(data[possition + 2]) << 8 | Byte.toUnsignedInt(data[possition + 3]);

            labels[i] = new AtomLabel();
            labels[i].atom = Byte.toUnsignedInt(data[possition + 5]) * 256 | Byte.toUnsignedInt(data[possition + 6]);
            labels[i].label = new byte[size - 7];

            for(int j = 0; j < size - 7; j++)
                labels[i].label[j] = data[possition + 7 + j];

            possition += size;
        }


        SGroup[] sgroups = null;

        if(sgroupCount > 0)
        {
            sgroups = new SGroup[sgroupCount];

            for(int i = 0; i < sgroupCount; i++)
            {
                sgroups[i] = new SGroup();

                sgroups[i].type = data[possition + 5];
                sgroups[i].subtype = data[possition + 6];
                sgroups[i].connectivity = data[possition + 7];

                int atomLength = Byte.toUnsignedInt(data[possition + 8]) * 256
                        | Byte.toUnsignedInt(data[possition + 9]);
                int bondLength = Byte.toUnsignedInt(data[possition + 10]) * 256
                        | Byte.toUnsignedInt(data[possition + 11]);

                possition += 12;

                sgroups[i].atoms = new int[atomLength];
                sgroups[i].bonds = new int[bondLength][];

                for(int a = 0; a < atomLength; a++)
                {
                    sgroups[i].atoms[a] = Byte.toUnsignedInt(data[possition + 0]) * 256
                            | Byte.toUnsignedInt(data[possition + 1]);
                    possition += 2;
                }

                for(int b = 0; b < bondLength; b++)
                {
                    sgroups[i].bonds[b] = new int[2];
                    sgroups[i].bonds[b][0] = Byte.toUnsignedInt(data[possition + 0]) * 256
                            | Byte.toUnsignedInt(data[possition + 1]);
                    sgroups[i].bonds[b][1] = Byte.toUnsignedInt(data[possition + 2]) * 256
                            | Byte.toUnsignedInt(data[possition + 3]);
                    possition += 4;
                }
            }
        }


        this.originalAtomCount = originalAtomCount;
        this.originalBondCount = originalBondCount;
        this.atomCount = atomCount;
        this.bondCount = bondCount;
        this.atomNumbers = atomNumbers;
        this.atomHydrogens = atomHydrogens;
        this.atomCharges = atomCharges;
        this.atomMasses = atomMasses;
        this.atomRadicalTypes = atomRadicalTypes;
        this.atomStereo = atomStereo;
        this.bondTypes = bondTypes;
        this.bondStereo = bondStereo;
        this.restH = restH;
        this.contains = contains;
        this.bondMatrix = bondMatrix;
        this.bondLists = new int[atomCount][];
        this.sgroups = sgroups;
        this.labels = labels;

        for(int i = 0; i < atomCount; i++)
        {
            ArrayList<Integer> tmpBondList = tmpBondLists.get(i);

            this.bondLists[i] = new int[tmpBondList.size()];

            for(int j = 0; j < tmpBondList.size(); j++)
                this.bondLists[i][j] = tmpBondList.get(j);
        }
    }


    public BinaryMolecule(byte[] data)
    {
        this(data, null, false, false, false, false, false, false, false, false, false);
    }


    @Override
    public final int getOriginalAtomCount()
    {
        return originalAtomCount;
    }


    @Override
    public final int getOriginalBondCount()
    {
        return originalBondCount;
    }


    @Override
    public final int getAtomCount()
    {
        return atomCount;
    }


    @Override
    public final int getBondCount()
    {
        return bondCount;
    }


    @Override
    public final boolean hasRestHydrogenFlags()
    {
        return restH != null;
    }


    @Override
    public final byte getAtomNumber(int atom)
    {
        return atomNumbers[atom];
    }


    @Override
    public AtomLabel getAtomLabel(int atom)
    {
        for(AtomLabel label : labels)
            if(label.atom == atom)
                return label;

        return null;
    }


    @Override
    public final byte getAtomHydrogenCount(int atom)
    {
        return atomHydrogens[atom];
    }


    @Override
    public final byte getAtomFormalCharge(int atom)
    {
        return atomCharges[atom];
    }


    @Override
    public final byte getAtomMass(int atom)
    {
        return atomMasses[atom];
    }


    @Override
    public final byte getAtomRadicalType(int atom)
    {
        return atomRadicalTypes[atom];
    }


    @Override
    public final byte getAtomStereo(int atom)
    {
        return atomStereo[atom];
    }


    @Override
    public final boolean getAtomRestHydrogenFlag(int atom)
    {
        return restH[atom];
    }


    @Override
    public final int getBond(int atom0, int atom1)
    {
        return bondMatrix[atom0 * atomCount + atom1];
    }


    @Override
    public final byte getBondType(int bond)
    {
        return bondTypes[bond];
    }


    @Override
    public final byte getBondStereo(int bond)
    {
        return bondStereo[bond];
    }


    @Override
    public final int getBondAtom(int bond, int i)
    {
        return contains[bond * 2 + i];
    }


    @Override
    public final boolean isAtomInBond(int atom, int bond)
    {
        return contains[bond * atomCount + 0] == atom || contains[bond * atomCount + 1] == atom;
    }


    @Override
    public final int[] getBondedAtoms(int atom)
    {
        return bondLists[atom];
    }


    public static boolean isExtended(byte[] data, boolean withRGroups, boolean withCharges, boolean withIsotopes,
            boolean withRadicals)
    {
        int possition = 0;

        int xAtomCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);
        int cAtomCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);
        int hAtomCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);
        int xBondCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);
        int specialCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);

        int heavyAtomCount = xAtomCount + cAtomCount;

        for(int i = 0; i < xAtomCount; i++)
            if(withRGroups && data[possition++] < 0)
                return true;


        int[] hBonds = new int[hAtomCount];

        for(int i = 0; i < hAtomCount; i++)
            hBonds[i] = 0;

        for(int i = 0; i < xBondCount; i++)
        {
            int b0 = Byte.toUnsignedInt(data[possition++]);
            int b1 = Byte.toUnsignedInt(data[possition++]);
            int b2 = Byte.toUnsignedInt(data[possition++]);
            possition++;

            int x = b0 | (b1 << 4 & 0xF00);
            int y = b2 | (b1 << 8 & 0xF00);

            if(x >= heavyAtomCount)
                hBonds[x - heavyAtomCount]++;

            if(y >= heavyAtomCount)
                hBonds[y - heavyAtomCount]++;
        }


        for(int i = 0; i < hAtomCount; i++)
        {
            int value = Byte.toUnsignedInt(data[possition++]) * 256 | Byte.toUnsignedInt(data[possition++]);

            if(value == 0 || (value & 0xFFF) >= heavyAtomCount)
                return true;

            hBonds[i]++;
        }


        for(int i = 0; i < hAtomCount; i++)
            if(hBonds[i] != 1)
                return true;

        if(!withCharges && !withIsotopes && !withRadicals)
            return false;

        for(int i = 0; i < specialCount; i++)
        {
            byte type = (byte) (Byte.toUnsignedInt(data[possition]) >> 4);
            int value = Byte.toUnsignedInt(data[possition++]) * 256 | Byte.toUnsignedInt(data[possition++]);
            possition++;

            int idx = value & 0xFFF;

            switch(type >> 4)
            {
                case SpecialRecordType.CHARGE:
                    if(withCharges && idx >= heavyAtomCount)
                        return true;
                    break;

                case SpecialRecordType.ISOTOPE:
                    if(withIsotopes && idx >= heavyAtomCount)
                        return true;
                    break;

                case SpecialRecordType.RADICAL:
                    if(withRadicals && idx >= heavyAtomCount)
                        return true;
                    break;
            }
        }

        return false;
    }


    public static boolean hasMultivalentHydrogen(byte[] data)
    {
        int xAtomCount = Byte.toUnsignedInt(data[0]) << 8 | Byte.toUnsignedInt(data[1]);
        int cAtomCount = Byte.toUnsignedInt(data[2]) << 8 | Byte.toUnsignedInt(data[3]);
        int xBondCount = Byte.toUnsignedInt(data[6]) << 8 | Byte.toUnsignedInt(data[7]);

        if(xBondCount == 0)
            return false;

        int offset = 10 + xAtomCount + (xBondCount - 1) * BOND_BLOCK_SIZE;

        int b0 = Byte.toUnsignedInt(data[offset + 0]);
        int b1 = Byte.toUnsignedInt(data[offset + 1]);
        int b2 = Byte.toUnsignedInt(data[offset + 2]);

        int x = b0 | (b1 << 4 & 0xF00);
        int y = b2 | (b1 << 8 & 0xF00);

        return x >= xAtomCount + cAtomCount || y >= xAtomCount + cAtomCount;
    }


    public static boolean hasPseudoAtom(byte[] data)
    {
        int possition = 0;

        int xAtomCount = Byte.toUnsignedInt(data[possition++]) << 8 | Byte.toUnsignedInt(data[possition++]);

        possition += 8;

        for(int i = 0; i < xAtomCount; i++)
            if((data[possition++]) < 0)
                return true;

        return false;
    }


    @Override
    public SGroup[] getSGroups()
    {
        return sgroups;
    }
}
