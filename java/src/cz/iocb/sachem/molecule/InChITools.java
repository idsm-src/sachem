package cz.iocb.sachem.molecule;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry.Conformation;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.ExtendedCisTrans;
import org.openscience.cdk.stereo.ExtendedTetrahedral;
import org.openscience.cdk.stereo.TetrahedralChirality;
import cz.iocb.sachem.molecule.Molecule.AtomType;



public class InChITools
{
    @SuppressWarnings("serial")
    public static class InChIException extends CDKException
    {
        public InChIException(String message)
        {
            super(message);
        }
    }


    public static class TautomericGroup
    {
        public int mobiles;
        public int charges;
        public List<Integer> atoms;
    }


    private static abstract class Field
    {
        private static final int elname = 0;
        private static final int el_number = 6;
        private static final int neighbor = 8;
        private static final int orig_at_number = 48;
        private static final int bond_stereo = 52;
        private static final int bond_type = 72;
        private static final int valence = 92;
        private static final int chem_bonds_valence = 93;
        private static final int num_H = 94;
        private static final int iso_atw_diff = 98;
        private static final int charge = 99;
        private static final int radical = 100;
        private static final int x = 112;
        private static final int y = 120;
        private static final int z = 128;
    }


    private static final int ZERO_ATW_DIFF = 127;
    private static final int MAXVAL = 20;
    private static final int RECORD_SIZE = 176;
    private static Isotopes isotopes;

    private final IAtomContainer molecule;

    @SuppressWarnings("rawtypes") private final List<IStereoElement> stereo = new ArrayList<IStereoElement>();
    private final Set<IAtom> stereoAtoms = new HashSet<IAtom>();
    private final Set<IBond> stereoBonds = new HashSet<IBond>();

    private final List<TautomericGroup> tautomericGroups = new ArrayList<TautomericGroup>();
    private final List<Integer> tautomericBonds = new ArrayList<Integer>();


    static
    {
        try
        {
            isotopes = Isotopes.getInstance();
        }
        catch(IOException e)
        {
        }
    }


    public InChITools(IAtomContainer molecule, boolean tautomers, boolean ignoreAnyBonds) throws CDKException
    {
        this.molecule = molecule;

        if(process(moleculeAsBytes(molecule, ignoreAnyBonds), tautomers) < 0)
            throw new InChIException("generation fails");

        molecule.setStereoElements(stereo);
    }


    private native int process(byte[] molecule, boolean tautomerism);


    void setStereoAtoms(short[] data)
    {
        for(int i = 0; i < data.length;)
        {
            final IAtom atom = molecule.getAtom(data[i++]);
            final int parity = data[i++];

            if(parity == -1 || parity == -2)
            {
                List<IBond> bonds = molecule.getConnectedBondsList(atom);
                IAtom[] peripherals = new IAtom[4];
                int position = 0;

                for(int b = 0; b < 2; b++)
                {
                    IAtom terminal = getTerminalAtom(atom, bonds.get(b));
                    List<IAtom> atoms = new ArrayList<IAtom>(2);

                    for(IBond bond : molecule.getConnectedBondsList(terminal))
                        if(bond.getOrder() != Order.DOUBLE)
                            atoms.add(bond.getOther(terminal));

                    if(atoms.size() == 1)
                        peripherals[position++] = terminal;

                    for(int diff = 0; diff < 3; diff++)
                        for(IAtom a : atoms)
                            if(a.getAtomicNumber() == AtomType.H && getMassDiff(a) == diff)
                                peripherals[position++] = a;

                    for(IAtom a : atoms)
                        if(!isInchiH(a))
                            peripherals[position++] = a;
                }

                stereo.add(new ExtendedTetrahedral(atom, peripherals, parity == -1 ?
                        ITetrahedralChirality.Stereo.CLOCKWISE : ITetrahedralChirality.Stereo.ANTI_CLOCKWISE));
            }
            else if(parity == 1 || parity == 2)
            {
                IAtom[] ligands = new IAtom[4];

                List<IAtom> neighbours = molecule.getConnectedAtomsList(atom);

                int position = 0;

                if(neighbours.size() == 3)
                    ligands[position++] = atom;

                for(int diff = 0; diff < 3; diff++)
                    for(IAtom a : neighbours)
                        if(a.getAtomicNumber() == AtomType.H && getMassDiff(a) == diff)
                            ligands[position++] = a;

                for(IAtom a : neighbours)
                    if(!isInchiH(a))
                        ligands[position++] = a;

                stereo.add(new TetrahedralChirality(atom, ligands,
                        parity == 1 ? Stereo.ANTI_CLOCKWISE : Stereo.CLOCKWISE));
            }

            stereoAtoms.add(atom);
        }
    }


    void setStereoBonds(short[] data)
    {
        for(int i = 0; i < data.length; i++)
        {
            final IBond bond = molecule.getBond(molecule.getAtom(data[i++]), molecule.getAtom(data[i++]));
            final int parity = data[i];

            if(parity == -1 || parity == -2)
            {
                IBond bonds[] = new IBond[2];

                for(int j = 0; j < 2; j++)
                {
                    IAtom terminal = getTerminalAtom(bond.getAtom(j), bond);
                    List<IBond> bondsList = molecule.getConnectedBondsList(terminal);

                    for(Iterator<IBond> it = bondsList.iterator(); it.hasNext();)
                        if(it.next().getOrder() == Order.DOUBLE)
                            it.remove();

                    IAtom other0 = bondsList.get(0).getOther(terminal);
                    IAtom other1 = bondsList.size() > 1 ? bondsList.get(1).getOther(terminal) : null;

                    if(other1 == null)
                        bonds[j] = bondsList.get(0);
                    else if(!isInchiH(other1))
                        bonds[j] = bondsList.get(1);
                    else if(!isInchiH(other0))
                        bonds[j] = bondsList.get(0);
                    else if(getMassDiff(other1) > getMassDiff(other0))
                        bonds[j] = bondsList.get(1);
                    else
                        bonds[j] = bondsList.get(0);
                }

                stereo.add(new ExtendedCisTrans(bond, bonds,
                        parity == -1 ? IStereoElement.TOGETHER : IStereoElement.OPPOSITE));
            }
            else if(parity == 1 || parity == 2)
            {
                IBond bonds[] = new IBond[2];

                for(int j = 0; j < 2; j++)
                {
                    IAtom atom = bond.getAtom(j);
                    List<IBond> bondsList = molecule.getConnectedBondsList(atom);
                    bondsList.remove(bond);

                    IAtom other0 = bondsList.get(0).getOther(atom);
                    IAtom other1 = bondsList.size() > 1 ? bondsList.get(1).getOther(atom) : null;

                    if(other1 == null)
                        bonds[j] = bondsList.get(0);
                    else if(!isInchiH(other1))
                        bonds[j] = bondsList.get(1);
                    else if(!isInchiH(other0))
                        bonds[j] = bondsList.get(0);
                    else if(getMassDiff(other1) > getMassDiff(other0))
                        bonds[j] = bondsList.get(1);
                    else
                        bonds[j] = bondsList.get(0);
                }

                stereo.add(new DoubleBondStereochemistry(bond, bonds,
                        parity == 2 ? Conformation.OPPOSITE : Conformation.TOGETHER));
            }

            stereoBonds.add(bond);
        }
    }


    void setAlternatingBonds(short[] data)
    {
        for(int i = 0; i < data.length; i++)
            tautomericBonds
                    .add(molecule.indexOf(molecule.getBond(molecule.getAtom(data[i++]), molecule.getAtom(data[i]))));
    }


    void setTautomericGroup(short[] data)
    {
        TautomericGroup group = new TautomericGroup();
        tautomericGroups.add(group);

        group.mobiles = data[0];
        group.charges = data[1];
        group.atoms = new ArrayList<Integer>(data.length - 2);


        for(int i = 2; i < data.length; i++)
            group.atoms.add((int) data[i]);
    }


    private final IAtom getTerminalAtom(IAtom atom, IBond bond)
    {
        while(true)
        {
            List<IBond> bonds = molecule.getConnectedBondsList(atom);
            bonds.remove(bond);

            if(bonds.size() != 1)
                return atom;

            bond = bonds.get(0);

            if(bond.getOrder() != Order.DOUBLE)
                return atom;

            atom = bond.getOther(atom);
        }
    }


    private static final byte[] moleculeAsBytes(IAtomContainer molecule, boolean ignoreAnyBonds) throws CDKException
    {
        if(molecule.getAtomCount() > Short.MAX_VALUE)
            throw new InChIException("too many atoms");

        int pseudoAtomId = Byte.MIN_VALUE;

        byte array[] = new byte[RECORD_SIZE * molecule.getAtomCount()];
        ByteBuffer buffer = ByteBuffer.wrap(array);
        buffer.order(ByteOrder.nativeOrder());

        for(int index = 0; index < molecule.getAtomCount(); index++)
        {
            int offset = index * RECORD_SIZE;
            IAtom atom = molecule.getAtom(index);
            List<IBond> bonds = molecule.getConnectedBondsList(atom);

            if(ignoreAnyBonds)
            {
                for(Iterator<IBond> it = bonds.iterator(); it.hasNext();)
                {
                    IBond bond = it.next();

                    if(bond instanceof QueryBond && ((QueryBond) bond).getExpression().type() == Expr.Type.TRUE)
                        it.remove();
                }
            }

            int multiplicity = molecule.getConnectedSingleElectronsCount(atom);

            int valence = 0;

            for(IBond bond : bonds)
                valence += getBondType(bond);

            if(atom.getImplicitHydrogenCount() + bonds.size() > MAXVAL)
                throw new InChIException("too large number of bonds");

            buffer.put(offset + Field.el_number, (byte) (int) atom.getAtomicNumber());
            buffer.putShort(offset + Field.orig_at_number, (short) index);
            buffer.put(offset + Field.valence, (byte) bonds.size());
            buffer.put(offset + Field.chem_bonds_valence, (byte) valence);
            buffer.put(offset + Field.num_H, (byte) (int) atom.getImplicitHydrogenCount());
            buffer.put(offset + Field.charge, (byte) (int) atom.getFormalCharge());
            buffer.put(offset + Field.radical, (byte) (multiplicity > 0 ? multiplicity + 1 : 0));

            if(!(atom instanceof IPseudoAtom))
            {
                byte[] name = atom.getSymbol().getBytes();

                for(int i = 0; i < name.length; i++)
                    buffer.put(offset + Field.elname + i, name[i]);

                int iso = getMassDiff(atom);
                buffer.put(offset + Field.iso_atw_diff, (byte) (iso == ZERO_ATW_DIFF ? 1 : iso > 0 ? iso + 1 : iso));
            }
            else if(pseudoAtomId < Byte.MAX_VALUE)
            {
                buffer.put(offset + Field.iso_atw_diff, (byte) (pseudoAtomId++));
            }
            else
            {
                throw new InChIException("too many pseudo atoms");
            }

            if(atom.getPoint3d() != null)
            {
                buffer.putDouble(offset + Field.x, atom.getPoint3d().x);
                buffer.putDouble(offset + Field.y, atom.getPoint3d().y);
                buffer.putDouble(offset + Field.z, atom.getPoint3d().z);
            }
            else if(atom.getPoint2d() != null)
            {
                buffer.putDouble(offset + Field.x, atom.getPoint2d().x);
                buffer.putDouble(offset + Field.y, atom.getPoint2d().y);
            }

            for(int i = 0; i < bonds.size(); i++)
            {
                IBond bond = bonds.get(i);

                buffer.putShort(offset + Field.neighbor + 2 * i, (short) molecule.indexOf(bond.getOther(atom)));
                buffer.put(offset + Field.bond_type + i, (byte) getBondType(bond));
                buffer.put(offset + Field.bond_stereo + i, (byte) getStereo(atom, bond));
            }
        }

        return array;
    }


    private static final int getMassDiff(IAtom atom)
    {
        if(atom.getMassNumber() == null)
            return 0;

        if(atom.getAtomicNumber() == 0)
            return 0;

        IIsotope majorIsotope = isotopes.getMajorIsotope(atom.getAtomicNumber());

        if(majorIsotope == null)
            return 0;

        return atom.getMassNumber() - majorIsotope.getMassNumber();
    }


    private static final int getBondType(IBond bond) throws CDKException
    {
        if(bond.getOrder() == Order.SINGLE)
            return 1;
        else if(bond.getOrder() == Order.DOUBLE)
            return 2;
        else if(bond.getOrder() == Order.TRIPLE)
            return 3;
        else
            throw new InChIException("unsupported bond type");
    }


    private static final int getStereo(IAtom atom, IBond bond)
    {
        int multiplier = bond.getAtom(0) == atom ? 1 : -1;

        switch(bond.getStereo())
        {
            case UP:
                return multiplier * 1;

            case UP_INVERTED:
                return multiplier * -1;

            case UP_OR_DOWN:
                return multiplier * 4;

            case UP_OR_DOWN_INVERTED:
                return multiplier * -4;

            case DOWN:
                return multiplier * 6;

            case DOWN_INVERTED:
                return multiplier * -6;

            case E_OR_Z:
                return 3;

            default:
                return 0;
        }
    }


    private static final boolean isInchiH(IAtom atom)
    {
        return atom.getAtomicNumber() == AtomType.H && getMassDiff(atom) < 3 && getMassDiff(atom) >= 0;
    }


    @SuppressWarnings("rawtypes")
    public List<IStereoElement> getStereoElements()
    {
        return stereo;
    }


    public Set<IAtom> getStereoAtoms()
    {
        return stereoAtoms;
    }


    public Set<IBond> getStereoBonds()
    {
        return stereoBonds;
    }


    public List<TautomericGroup> getTautomericGroups()
    {
        return tautomericGroups;
    }


    public List<Integer> getTautomericBonds()
    {
        return tautomericBonds;
    }
}
