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
    private static final boolean lazyInitialization = true;

    private static final String idFieldName = "id";
    private static final String subfpFieldName = "subfp";
    private static final String simfpFieldName = "simfp";
    private static final String simSizeFieldName = "simsz";

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
                builder.add(new TermQuery(new Term(subfpFieldName, tokenizer.bitAsString(bit))),
                        BooleanClause.Occur.MUST);
        }
        else
        {
            for(int bit : fp)
                builder.add(IntPoint.newExactQuery(subfpFieldName, bit), BooleanClause.Occur.MUST);
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
                builder.add(new TermQuery(new Term(simfpFieldName, tokenizer.bitAsString(bit))),
                        BooleanClause.Occur.SHOULD);
        }
        else
        {
            for(int bit : fp)
                builder.add(IntPoint.newExactQuery(simfpFieldName, bit), BooleanClause.Occur.SHOULD);
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


            if(lazyInitialization)
            {
                if(useIdTable)
                {
                    for(int i = 0; i < idTable.length; i++)
                        idTable[i] = Integer.MIN_VALUE;
                }

                if(useSizeTable)
                {
                    for(int i = 0; i < sizeTable.length; i++)
                        sizeTable[i] = Integer.MIN_VALUE;
                }
            }
            else
            {
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
                            StoredField field = (StoredField) doc.getField(simSizeFieldName);
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
    }


    protected final int getMoleculeId(int id) throws IOException
    {
        if(useIdTable)
        {
            if(lazyInitialization && idTable[id] == Integer.MIN_VALUE)
            {
                Document doc = searcher.doc(id);
                IndexableField field = doc.getField(idFieldName);
                idTable[id] = field.numericValue().intValue();
            }

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


    protected final int getMoleculeSimFpSize(int id) throws IOException
    {
        if(useSizeTable)
        {
            if(lazyInitialization && sizeTable[id] == Integer.MIN_VALUE)
            {
                Document document = searcher.doc(id);
                StoredField field = (StoredField) document.getField(simSizeFieldName);
                sizeTable[id] = field.numericValue().intValue();
            }

            return sizeTable[id];
        }
        else
        {
            Document document = searcher.doc(id);
            StoredField field = (StoredField) document.getField(simSizeFieldName);
            int fpSize = field.numericValue().intValue();

            return fpSize;
        }
    }
}
