package cz.iocb.sachem.search;

import java.nio.charset.StandardCharsets;
import java.util.concurrent.atomic.AtomicInteger;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.sachem.isomorphism.IsomorphismSort;
import cz.iocb.sachem.shared.MoleculeCreator;



public class LucyLoader
{
    public static class LucyData
    {
        public byte[] molecule;
        public String exception;
    }


    public static LucyData[] getIndexData(final byte[][] molfilesArray) throws InterruptedException
    {
        final int cores = Runtime.getRuntime().availableProcessors();
        final AtomicInteger idx = new AtomicInteger(0);
        final LucyData[] result = new LucyData[molfilesArray.length];

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
                        LucyData item = new LucyData();

                        try
                        {
                            String molfile = new String(molfilesArray[i], StandardCharsets.ISO_8859_1);

                            IAtomContainer readMolecule = MoleculeCreator.getMoleculeFromMolfile(molfile);
                            MoleculeCreator.configureMolecule(readMolecule);

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
}
