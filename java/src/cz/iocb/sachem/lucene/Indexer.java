package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import org.apache.lucene.document.BinaryDocValuesField;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.document.NumericDocValuesField;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.SerialMergeScheduler;
import org.apache.lucene.search.similarities.BooleanSimilarity;
import org.apache.lucene.store.FSDirectory;
import org.apache.lucene.util.BytesRef;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.sachem.fingerprint.IOCBFingerprint;
import cz.iocb.sachem.molecule.AromaticityMode;
import cz.iocb.sachem.molecule.BinaryMolecule;
import cz.iocb.sachem.molecule.BinaryMoleculeBuilder;
import cz.iocb.sachem.molecule.MoleculeCreator;



public class Indexer
{
    private static class IndexItem
    {
        private int id;
        private byte[] molecule;

        private IndexItem(int id, byte[] molecule)
        {
            this.id = id;
            this.molecule = molecule;
        }
    }


    private FSDirectory folder;
    private IndexWriter indexer;

    private Thread[] documentThreads;
    private Thread indexThread;
    private ArrayBlockingQueue<IndexItem> moleculeQueue;
    private ArrayBlockingQueue<Document> documentQueue;


    private Throwable exception;


    public void begin(String path, int segments, int bufferedDocs, double bufferSize) throws IOException
    {
        folder = FSDirectory.open(Paths.get(path));

        if(bufferedDocs <= 0)
            bufferedDocs = IndexWriterConfig.DISABLE_AUTO_FLUSH;

        if(bufferSize <= 0)
            bufferSize = IndexWriterConfig.DISABLE_AUTO_FLUSH;

        try
        {
            BalancedMergePolicy policy = new BalancedMergePolicy(segments);

            IndexWriterConfig config = new IndexWriterConfig(null);

            config.setMaxBufferedDocs(bufferedDocs);
            config.setRAMBufferSizeMB(bufferSize);

            config.setMergePolicy(policy);
            config.setSimilarity(new BooleanSimilarity());
            config.setMergeScheduler(new SerialMergeScheduler());

            final int cores = Runtime.getRuntime().availableProcessors();
            indexer = new IndexWriter(folder, config);

            moleculeQueue = new ArrayBlockingQueue<IndexItem>(128 * cores);
            documentQueue = new ArrayBlockingQueue<Document>(2 * bufferedDocs);

            documentThreads = new Thread[cores];

            for(int i = 0; i < cores; i++)
            {
                Thread thread = new Thread()
                {
                    @Override
                    public void run()
                    {
                        while(true)
                        {
                            try
                            {
                                IndexItem item = moleculeQueue.take();

                                if(item.molecule == null)
                                    return;

                                if(exception != null)
                                    continue;

                                Document document = createDocument(item.id, item.molecule);

                                synchronized(indexer)
                                {
                                    indexer.addDocument(document);
                                }
                            }
                            catch(Throwable e)
                            {
                                synchronized(Indexer.this)
                                {
                                    if(exception != null)
                                        exception = e;
                                }
                            }
                        }
                    }
                };

                thread.start();

                documentThreads[i] = thread;
            }


            Thread thread = new Thread()
            {
                @Override
                public void run()
                {
                    try
                    {
                        while(true)
                        {
                            Document document = documentQueue.take();

                            if(document.getFields().isEmpty())
                                break;

                            indexer.addDocument(document);
                        }
                    }
                    catch(Throwable e)
                    {
                        synchronized(Indexer.this)
                        {
                            if(exception != null)
                                exception = e;
                        }
                    }
                }
            };

            thread.start();

            indexThread = thread;
        }
        catch(Throwable e)
        {
            stopThreads();

            folder.close();
            throw e;
        }
    }


    public String add(int id, byte[] molfile) throws IOException, InterruptedException
    {
        if(exception != null)
            throw new IOException(exception);

        try
        {
            IAtomContainer container = MoleculeCreator.getMoleculeFromMolfile(new String(molfile),
                    AromaticityMode.AUTO);
            BinaryMoleculeBuilder builder = new BinaryMoleculeBuilder(container);

            byte[] binary = builder.asBytes(true);
            moleculeQueue.put(new IndexItem(id, binary));
        }
        catch(Exception e)
        {
            return e.getMessage();
        }

        return null;
    }


    public void delete(int id) throws IOException
    {
        indexer.deleteDocuments(IntPoint.newExactQuery(Settings.idFieldName, id));
    }


    public void optimize() throws IOException
    {
        indexer.forceMergeDeletes();
    }


    public void commit() throws IOException
    {
        stopThreads();

        if(exception != null)
            throw new IOException(exception);

        indexer.commit();
        indexer.close();
        folder.close();
    }


    public void rollback() throws IOException
    {
        stopThreads();

        indexer.rollback();
        indexer.close();
        folder.close();
    }


    private void stopThreads()
    {
        try
        {
            if(documentThreads != null)
            {
                for(int i = 0; i < documentThreads.length; i++)
                    moleculeQueue.put(new IndexItem(Integer.MIN_VALUE, null));

                for(int i = 0; i < documentThreads.length; i++)
                    if(documentThreads[i] != null)
                        documentThreads[i].join();
            }

            if(indexThread != null)
            {
                documentQueue.put(new Document());
                indexThread.join();
            }
        }
        catch(InterruptedException e)
        {
        }
    }


    private static Document createDocument(int id, byte[] binary)
    {
        Document document = new Document();
        document.add(new IntPoint(Settings.idFieldName, id));
        document.add(new NumericDocValuesField(Settings.idFieldName, id));
        document.add(new StoredField(Settings.idFieldName, id));

        BinaryMolecule molecule = new BinaryMolecule(binary);


        /* substructure index */
        Set<Integer> subFp = IOCBFingerprint.getSubstructureFingerprint(molecule);

        document.add(new StoredField(Settings.substructureFieldName, binary));
        document.add(new BinaryDocValuesField(Settings.substructureFieldName, new BytesRef(binary)));
        document.add(new IntPoint(Settings.substructureFieldName, subFp.size()));
        document.add(new TextField(Settings.substructureFieldName, new FingerprintTokenStream(subFp)));


        /* similarity index */
        List<List<Integer>> simFp = IOCBFingerprint.getSimilarityFingerprint(molecule, Settings.maximumSimilarityDepth);

        byte[] array = new byte[simFp.stream().map(i -> i.size() + 1).reduce(0, Integer::sum) * Integer.BYTES];

        for(int pos = 0, i = 0; i < simFp.size(); i++)
        {
            for(int b = 0; b < Integer.BYTES; b++)
                array[pos++] = (byte) (simFp.get(i).size() >> (8 * b));

            for(int bit : simFp.get(i))
                for(int b = 0; b < Integer.BYTES; b++)
                    array[pos++] = (byte) (bit >> (8 * b));
        }


        Set<Integer> bits = new HashSet<Integer>();

        for(List<Integer> seg : simFp)
            bits.addAll(seg);


        document.add(new TextField(Settings.similarityFieldName, new FingerprintTokenStream(bits)));
        document.add(new StoredField(Settings.similarityFieldName, array));
        document.add(new BinaryDocValuesField(Settings.similarityFieldName, new BytesRef(array)));

        for(int size = 0, i = 0; i < simFp.size(); i++)
        {
            size += simFp.get(i).size();
            document.add(new IntPoint(Settings.similarityFieldName, size));
            size += SimilarStructureQuery.iterationSizeOffset;
        }

        return document;
    }
}
