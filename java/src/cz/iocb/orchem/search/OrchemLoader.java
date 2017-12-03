package cz.iocb.orchem.search;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.orchem.fingerprint.OrchemExtendedFingerprinter;
import cz.iocb.orchem.isomorphism.IsomorphismSort;
import cz.iocb.orchem.shared.MoleculeCounts;
import cz.iocb.orchem.shared.MoleculeCreator;



public class OrchemLoader
{
    public static class OrchemData
    {
        public long[] fp;
        public short[] counts;
        public byte[] atoms;
        public byte[] bonds;
    }


    private static final ThreadLocal<OrchemExtendedFingerprinter> fingerPrinter = new ThreadLocal<OrchemExtendedFingerprinter>()
    {
        @Override
        protected OrchemExtendedFingerprinter initialValue()
        {
            return new OrchemExtendedFingerprinter();
        }
    };


    public static OrchemData getIndexData(byte[] molfileArray) throws CDKException, IOException
    {
        String molfile = new String(molfileArray, StandardCharsets.ISO_8859_1);
        OrchemData item = new OrchemData();

        IAtomContainer readMolecule = MoleculeCreator.getMoleculeFromMolfile(molfile);
        MoleculeCreator.configureMolecule(readMolecule);

        // calculate similarity fingerprint
        BitSet fp = fingerPrinter.get().getBitFingerprint(readMolecule).asBitSet();
        item.fp = fp.toLongArray();

        // calculate molecule couts
        MoleculeCounts counts = new MoleculeCounts(readMolecule);
        item.counts = new short[10];
        item.counts[0] = counts.molTripleBondCount;
        item.counts[1] = counts.molSCount;
        item.counts[2] = counts.molOCount;
        item.counts[3] = counts.molNCount;
        item.counts[4] = counts.molFCount;
        item.counts[5] = counts.molClCount;
        item.counts[6] = counts.molBrCount;
        item.counts[7] = counts.molICount;
        item.counts[8] = counts.molCCount;
        item.counts[9] = counts.molPCount;

        // calculate molecule binary representation
        readMolecule.setAtoms(IsomorphismSort.atomsByFrequency(readMolecule));
        OrchemMoleculeBuilder builder = new OrchemMoleculeBuilder(readMolecule);
        item.atoms = builder.atomsAsBytes();
        item.bonds = builder.bondsAsBytes();

        return item;
    }
}
