package cz.iocb.sachem.molecule;

public abstract class Molecule
{
    public static abstract class AtomType
    {
        public static final byte H = 1;
        public static final byte C = 6;
        public static final byte N = 7;
        public static final byte O = 8;
        public static final byte F = 9;
        public static final byte P = 15;
        public static final byte S = 16;
        public static final byte Cl = 17;
        public static final byte Br = 35;
        public static final byte I = 53;

        public static final byte R = -'R';
        public static final byte Q = -'Q';
        public static final byte M = -'M';
        public static final byte X = -'X';

        public static final byte A = -'A';
        public static final byte G = -'G';

        public static final byte POSITRONIUM = -'p';
        public static final byte ELECTRON = -'e';
        public static final byte PHOTON = -'h';
        public static final byte MUONIUM = -'m';

        public static final byte ENZYME = -'z';
        public static final byte ACP = -'a';

        public static final byte EMPTY = -' ';
        public static final byte UNKNOWN = -'?';
    }


    public static abstract class BondType
    {
        public static final byte NONE = 0;
        public static final byte SINGLE = 1;
        public static final byte DOUBLE = 2;
        public static final byte TRIPLE = 3;
        public static final byte QUADRUPLE = 4;
        public static final byte QUINTUPLE = 5;
        public static final byte SEXTUPLE = 6;
        public static final byte AROMATIC = 11;
        public static final byte SINGLE_OR_DOUBLE = 12;
        public static final byte SINGLE_OR_AROMATIC = 13;
        public static final byte DOUBLE_OR_AROMATIC = 14;
        public static final byte ANY = 15;
    }


    public static abstract class TetrahedralStereo
    {
        public static final byte NONE = 0;
        public static final byte CLOCKWISE = 1;
        public static final byte ANTI_CLOCKWISE = 2;
        public static final byte UNDEFINED = 3;
    }


    public static abstract class BondStereo
    {
        public static final byte NONE = 0;
        public static final byte OPPOSITE = 1;
        public static final byte TOGETHER = 2;
        public static final byte UNDEFINED = 3;
    }


    public static final int MAX_ATOM_IDX = Short.MAX_VALUE;
}
