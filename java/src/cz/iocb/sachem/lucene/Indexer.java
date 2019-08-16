package cz.iocb.sachem.lucene;

import static cz.iocb.sachem.lucene.Settings.idFieldName;
import static cz.iocb.sachem.lucene.Settings.indexType;
import static cz.iocb.sachem.lucene.Settings.simSizeFieldName;
import static cz.iocb.sachem.lucene.Settings.simfpFieldName;
import static cz.iocb.sachem.lucene.Settings.subfpFieldName;
import java.io.IOException;
import java.nio.file.Paths;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.SerialMergeScheduler;
import org.apache.lucene.index.TieredMergePolicy;
import org.apache.lucene.search.similarities.BooleanSimilarity;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;



public class Indexer
{
    private FSDirectory folder;
    private IndexWriter indexer;


    public void begin(String path) throws IOException
    {
        folder = FSDirectory.open(Paths.get(path));

        try
        {
            TieredMergePolicy policy = new TieredMergePolicy();
            policy.setMaxMergeAtOnceExplicit(Integer.MAX_VALUE);
            policy.setMaxMergedSegmentMB(1024 * 1024);

            FingerprintAnalyzer analyzer = new FingerprintAnalyzer();
            IndexWriterConfig config = new IndexWriterConfig(analyzer);
            config.setMergePolicy(policy);
            config.setSimilarity(new BooleanSimilarity());
            config.setMergeScheduler(new SerialMergeScheduler());

            indexer = new IndexWriter(folder, config);
        }
        catch(Throwable e)
        {
            folder.close();
            throw e;
        }
    }


    public void add(int id, int[] subfp, int[] simfp) throws IOException
    {
        Document document = new Document();
        document.add(new IntPoint(idFieldName, id));
        document.add(new StoredField(idFieldName, id));
        document.add(new StoredField(simSizeFieldName, simfp.length));

        if(indexType == IndexType.TEXT)
        {
            document.add(new TextField(subfpFieldName, new FingerprintReader(subfp)));
            document.add(new TextField(simfpFieldName, new FingerprintReader(simfp)));
        }
        else
        {
            for(int bit : subfp)
                document.add(new IntPoint(subfpFieldName, bit));

            for(int bit : simfp)
                document.add(new IntPoint(simfpFieldName, bit));
        }

        indexer.addDocument(document);
    }


    public void addIndex(String path) throws IOException
    {
        Directory subFolder = FSDirectory.open(Paths.get(path));
        indexer.addIndexes(subFolder);
    }


    public void delete(int id) throws IOException
    {
        indexer.deleteDocuments(IntPoint.newExactQuery(idFieldName, id));
    }


    public void optimize() throws IOException
    {
        indexer.forceMerge(1);
    }


    public void commit() throws IOException
    {
        indexer.commit();
        indexer.close();
        folder.close();
    }


    public void rollback() throws IOException
    {
        indexer.rollback();
        indexer.close();
        folder.close();
    }
}
