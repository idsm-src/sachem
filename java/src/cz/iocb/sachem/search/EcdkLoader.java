package cz.iocb.sachem.search;

import java.nio.charset.StandardCharsets;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicInteger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.isomorphism.matchers.CTFileQueryBond;
import cz.iocb.sachem.fingerprint.ExtendedFingerprinter;
import cz.iocb.sachem.isomorphism.IsomorphismSort;
import cz.iocb.sachem.shared.MoleculeCreator;



public class EcdkLoader
{
    public static class EcdkData
    {
        public long[] fp;
        //public short[] counts;
        public byte[] molecule;
        public String exception;
    }


    private static final ThreadLocal<ExtendedFingerprinter> fingerPrinter = new ThreadLocal<ExtendedFingerprinter>()
    {
        @Override
        protected ExtendedFingerprinter initialValue()
        {
            return new ExtendedFingerprinter();
        }
    };


    public static EcdkData[] getIndexData(final byte[][] molfilesArray) throws InterruptedException
    {
        final int cores = Runtime.getRuntime().availableProcessors();
        final AtomicInteger idx = new AtomicInteger(0);
        final EcdkData[] result = new EcdkData[molfilesArray.length];

        Thread[] thread = new Thread[cores];

        for(int c = 0; c < cores; c++)
        {
            thread[c] = new Thread()
            {
                @Override
                public void run()
                {
                    for(int i = idx.getAndIncrement(); i < molfilesArray.length; i = idx.getAndIncrement())
                    {
                        EcdkData item = new EcdkData();

                        try
                        {
                            String molfile = new String(molfilesArray[i], StandardCharsets.ISO_8859_1);

                            IAtomContainer readMolecule = MoleculeCreator.getMoleculeFromMolfile(molfile);
                            MoleculeCreator.configureMolecule(readMolecule);

                            // calculate similarity fingerprint
                            if(!canBeIndexed(readMolecule))
                                throw new CDKException("query molecule cannot be indexed");

                            // calculate similarity fingerprint
                            BitSet fp = fingerPrinter.get().getFingerprint(readMolecule);
                            item.fp = fp.toLongArray();

                            // calculate molecule couts
                            /*
                            MoleculeCounts counts = new MoleculeCounts(readMolecule, true);
                            item.counts = new short[13];
                            item.counts[0] = counts.molSingleBondCount;
                            item.counts[1] = counts.molDoubleBondCount;
                            item.counts[2] = counts.molTripleBondCount;
                            item.counts[3] = counts.molAromaticBondCount;
                            item.counts[4] = counts.molSCount;
                            item.counts[5] = counts.molOCount;
                            item.counts[6] = counts.molNCount;
                            item.counts[7] = counts.molFCount;
                            item.counts[8] = counts.molClCount;
                            item.counts[9] = counts.molBrCount;
                            item.counts[10] = counts.molICount;
                            item.counts[11] = counts.molCCount;
                            item.counts[12] = counts.molPCount;
                            */

                            // calculate molecule binary representation
                            readMolecule.setAtoms(IsomorphismSort.atomsByFrequency(readMolecule));
                            SachemMoleculeBuilder builder = new SachemMoleculeBuilder(readMolecule);
                            item.molecule = builder.asBytes(true);
                        }
                        catch(Throwable e)
                        {
                            item.exception = "error: " + e.getClass().getCanonicalName() + ": " + e.getMessage();
                        }

                        result[i] = item;
                    }

                }
            };
        }

        for(int c = 0; c < cores; c++)
            thread[c].start();

        for(int c = 0; c < cores; c++)
            thread[c].join();

        return result;
    }


    private static boolean canBeIndexed(IAtomContainer molecule)
    {
        for(IAtom atom : molecule.atoms())
            if(atom instanceof IPseudoAtom)
                return false;

        for(IBond bond : molecule.bonds())
            if(bond instanceof CTFileQueryBond && !bond.isAromatic())
                return false;

        return true;
    }
}
