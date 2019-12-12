package cz.iocb.sachem.molecule;

public final class BinaryMolecule extends Molecule
{
    public static abstract class SpecialRecordType
    {
        public static final byte CHARGE = 0;
        public static final byte ISOTOPE = 1;
        public static final byte TETRAHEDRAL_STEREO = 2;
        public static final byte BOND_STEREO = 3;
    }
}
