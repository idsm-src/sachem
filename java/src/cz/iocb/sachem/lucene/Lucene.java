package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.BitSet;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.IndexableField;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.index.Term;
import org.apache.lucene.index.TieredMergePolicy;
import org.apache.lucene.search.BooleanClause;
import org.apache.lucene.search.BooleanQuery;
import org.apache.lucene.search.BooleanQuery.Builder;
import org.apache.lucene.search.ConstantScoreQuery;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.SimpleCollector;
import org.apache.lucene.search.TermQuery;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;



public class Lucene
{
    private static enum IndexType
    {
        TEXT, POINTS
    };


    private class BitSetCollector extends SimpleCollector
    {
        private final BitSet bitset;
        int docBase;

        BitSetCollector()
        {
            this.bitset = new BitSet(maxMoleculeId);
        }

        @Override
        protected void doSetNextReader(LeafReaderContext context) throws IOException
        {
            docBase = context.docBase;
        }

        @Override
        public void collect(int docId) throws IOException
        {
            Document doc = searcher.doc(docBase + docId);
            IndexableField id = doc.getField(idFieldName);
            int moleculeID = id.numericValue().intValue();

            bitset.set(moleculeID);
        }

        @Override
        public boolean needsScores()
        {
            return false;
        }

        public BitSet getBitSet()
        {
            return bitset;
        }
    }


    private static final IndexType indexType = IndexType.POINTS;
    private static final String idFieldName = "id";
    private static final String fpFieldName = "fp";

    private int maxMoleculeId;
    private Directory folder;
    private IndexSearcher searcher;
    private IndexWriter indexer;
    private final FingerprintTokenizer tokenizer = new FingerprintTokenizer();


    void setFolder(String path, int maxId) throws IOException
    {
        maxMoleculeId = maxId;
        folder = FSDirectory.open(Paths.get(path));
        searcher = null;
    }


    void begin() throws IOException
    {
        TieredMergePolicy policy = new TieredMergePolicy();
        policy.setMaxMergeAtOnceExplicit(Integer.MAX_VALUE);
        policy.setMaxMergedSegmentMB(1024 * 1024);

        FingerprintAnalyzer analyzer = new FingerprintAnalyzer();
        IndexWriterConfig config = new IndexWriterConfig(analyzer);
        config.setMergePolicy(policy);
        indexer = new IndexWriter(folder, config);
    }


    void add(int id, int[] fp) throws IOException
    {
        Document document = new Document();
        document.add(new IntPoint(idFieldName, id));
        document.add(new StoredField(idFieldName, id));

        if(indexType == IndexType.TEXT)
        {
            document.add(new TextField(fpFieldName, new FingerprintReader(fp)));
        }
        else
        {
            for(int bit : fp)
                document.add(new IntPoint(fpFieldName, bit));
        }

        indexer.addDocument(document);
    }


    void addIndex(String path) throws IOException
    {
        Directory subFolder = FSDirectory.open(Paths.get(path));
        indexer.addIndexes(subFolder);
    }


    void delete(int id) throws IOException
    {
        indexer.deleteDocuments(IntPoint.newExactQuery(idFieldName, id));
    }


    void optimize() throws IOException
    {
        indexer.forceMerge(1);
    }


    void commit() throws IOException
    {
        indexer.commit();
        indexer.close();
    }


    void rollback() throws IOException
    {
        indexer.rollback();
        indexer.close();
    }


    long[] search(int fp[], int maxResults) throws IOException
    {
        if(searcher == null)
        {
            IndexReader reader = DirectoryReader.open(folder);
            searcher = new IndexSearcher(reader);
        }


        Builder builder = new BooleanQuery.Builder();

        if(indexType == IndexType.TEXT)
        {
            for(int bit : fp)
                builder.add(new TermQuery(new Term(fpFieldName, tokenizer.bitAsString(bit))), BooleanClause.Occur.MUST);
        }
        else
        {
            for(int bit : fp)
                builder.add(IntPoint.newExactQuery(fpFieldName, bit), BooleanClause.Occur.MUST);
        }

        BitSetCollector collector = new BitSetCollector();
        searcher.search(new ConstantScoreQuery(builder.build()), collector);

        return collector.getBitSet().toLongArray();
    }
}
