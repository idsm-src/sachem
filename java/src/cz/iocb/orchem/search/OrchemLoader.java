package cz.iocb.orchem.search;

import java.nio.charset.StandardCharsets;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicInteger;
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
        public byte[] molecule;
        public String exception;
    }


    private static final ThreadLocal<OrchemExtendedFingerprinter> fingerPrinter = new ThreadLocal<OrchemExtendedFingerprinter>()
    {
        @Override
        protected OrchemExtendedFingerprinter initialValue()
        {
            return new OrchemExtendedFingerprinter();
        }
    };


    public static OrchemData[] getIndexData(final byte[][] molfilesArray) throws InterruptedException
    {
        final int cores = Runtime.getRuntime().availableProcessors();
        final AtomicInteger idx = new AtomicInteger(0);
        final OrchemData[] result = new OrchemData[molfilesArray.length];

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
                        OrchemData item = new OrchemData();

                        try
                        {
                            String molfile = new String(molfilesArray[i], StandardCharsets.ISO_8859_1);

                            IAtomContainer readMolecule = MoleculeCreator.getMoleculeFromMolfile(molfile);
                            MoleculeCreator.configureMolecule(readMolecule);

                            // calculate similarity fingerprint
                            BitSet fp = fingerPrinter.get().getFingerprint(readMolecule);
                            item.fp = fp.toLongArray();

                            // calculate molecule couts
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

                            // calculate molecule binary representation
                            readMolecule.setAtoms(IsomorphismSort.atomsByFrequency(readMolecule));
                            OrchemMoleculeBuilder builder = new OrchemMoleculeBuilder(readMolecule);
                            item.molecule = builder.asBytes(true);
                        }
                        catch(Throwable e)
                        {
                            item.exception = e.getClass().getCanonicalName() + ": " + e.getMessage();
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
}
