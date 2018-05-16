package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.nio.file.Paths;
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
import org.apache.lucene.search.MatchAllDocsQuery;
import org.apache.lucene.search.SimpleCollector;
import org.apache.lucene.search.TermQuery;
import org.apache.lucene.search.similarities.BooleanSimilarity;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;



public class Lucene
{
    private static enum IndexType
    {
        TEXT, POINTS
    };


    private static final IndexType indexType = IndexType.TEXT;
    private static final boolean useIdTable = true;
    private static final boolean useSizeTable = true;

    private static final String idFieldName = "id";
    private static final String fpFieldName = "fp";
    private static final String sizeFieldName = "sz";

    private Directory folder;
    private IndexSearcher searcher;
    private IndexWriter indexer;
    private int[] idTable;
    private int[] sizeTable;
    private final FingerprintTokenizer tokenizer = new FingerprintTokenizer();


    public void setFolder(String path) throws IOException
    {
        folder = FSDirectory.open(Paths.get(path));
        searcher = null;

        if(useIdTable)
            idTable = null;

        if(useSizeTable)
            sizeTable = null;
    }


    public void begin() throws IOException
    {
        TieredMergePolicy policy = new TieredMergePolicy();
        policy.setMaxMergeAtOnceExplicit(Integer.MAX_VALUE);
        policy.setMaxMergedSegmentMB(1024 * 1024);

        FingerprintAnalyzer analyzer = new FingerprintAnalyzer();
        IndexWriterConfig config = new IndexWriterConfig(analyzer);
        config.setMergePolicy(policy);
        config.setSimilarity(new BooleanSimilarity());
        indexer = new IndexWriter(folder, config);
    }


    public void add(int id, int[] fp) throws IOException
    {
        Document document = new Document();
        document.add(new IntPoint(idFieldName, id));
        document.add(new StoredField(idFieldName, id));
        document.add(new StoredField(sizeFieldName, fp.length));

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
    }


    public void rollback() throws IOException
    {
        indexer.rollback();
        indexer.close();
    }


    public long[] subsearch(int fp[], int maxMoleculeId) throws IOException
    {
        if(searcher == null)
            initSearcher();


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

        BitSetCollector collector = new BitSetCollector(this, maxMoleculeId);
        searcher.search(new ConstantScoreQuery(builder.build()), collector);

        return collector.getBitSet().toLongArray();
    }


    ScoreHit[] simsearch(int fp[], int top, float cutoff) throws IOException
    {
        if(searcher == null)
            initSearcher();


        Builder builder = new BooleanQuery.Builder();

        if(indexType == IndexType.TEXT)
        {
            for(int bit : fp)
                builder.add(new TermQuery(new Term(fpFieldName, tokenizer.bitAsString(bit))),
                        BooleanClause.Occur.SHOULD);
        }
        else
        {
            for(int bit : fp)
                builder.add(IntPoint.newExactQuery(fpFieldName, bit), BooleanClause.Occur.SHOULD);
        }


        SimilarDocCollector collector = new SimilarDocCollector(this, fp.length, top, cutoff);
        searcher.search(builder.build(), collector);

        return collector.getDocs();
    }


    protected void initSearcher() throws IOException
    {
        IndexReader reader = DirectoryReader.open(folder);
        searcher = new IndexSearcher(reader);
        searcher.setSimilarity(new BooleanSimilarity());

        @SuppressWarnings("unused")
        boolean useCache = useIdTable || useSizeTable;

        if(useCache)
        {
            if(useIdTable)
                idTable = new int[reader.maxDoc()];

            if(useSizeTable)
                sizeTable = new int[reader.maxDoc()];

            searcher.search(new MatchAllDocsQuery(), new SimpleCollector()
            {
                private int docBase;

                @Override
                protected void doSetNextReader(LeafReaderContext context) throws IOException
                {
                    docBase = context.docBase;
                }

                @Override
                public void collect(int docId) throws IOException
                {
                    Document doc = searcher.doc(docBase + docId);

                    if(useIdTable)
                    {
                        IndexableField id = doc.getField(idFieldName);
                        int moleculeId = id.numericValue().intValue();
                        idTable[docBase + docId] = moleculeId;
                    }

                    if(useSizeTable)
                    {
                        StoredField field = (StoredField) doc.getField(sizeFieldName);
                        int fpSize = field.numericValue().intValue();
                        sizeTable[docBase + docId] = fpSize;
                    }
                }

                @Override
                public boolean needsScores()
                {
                    return false;
                }
            });
        }
    }


    protected final int getMoleculeId(int id) throws IOException
    {
        if(useIdTable)
        {
            return idTable[id];
        }
        else
        {
            Document doc = searcher.doc(id);
            IndexableField field = doc.getField(idFieldName);
            int moleculeId = field.numericValue().intValue();

            return moleculeId;
        }
    }


    protected final int getMoleculeFpSize(int id) throws IOException
    {
        if(useSizeTable)
        {
            return sizeTable[id];
        }
        else
        {
            Document document = searcher.doc(id);
            StoredField field = (StoredField) document.getField(sizeFieldName);
            int fpSize = field.numericValue().intValue();

            return fpSize;
        }
    }
}
