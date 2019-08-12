package cz.iocb.sachem.lucene;

import static cz.iocb.sachem.lucene.Settings.idFieldName;
import static cz.iocb.sachem.lucene.Settings.indexType;
import static cz.iocb.sachem.lucene.Settings.lazyInitialization;
import static cz.iocb.sachem.lucene.Settings.simSizeFieldName;
import static cz.iocb.sachem.lucene.Settings.simfpFieldName;
import static cz.iocb.sachem.lucene.Settings.subfpFieldName;
import static cz.iocb.sachem.lucene.Settings.useIdTable;
import static cz.iocb.sachem.lucene.Settings.useSizeTable;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.IndexableField;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.index.Term;
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
    private Directory folder;
    private IndexSearcher searcher;
    private int[] idTable;
    private int[] sizeTable;
    private final FingerprintTokenizer tokenizer = new FingerprintTokenizer();


    public synchronized void setFolder(String pathName) throws IOException
    {
        close();

        Path path = Paths.get(pathName);
        folder = FSDirectory.open(path);
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


    public synchronized long[] subsearch(int fp[], int maxMoleculeId) throws IOException
    {
        if(fp.length == 0)
        {
            BitSetCollector collector = new BitSetCollector(this, maxMoleculeId);
            searcher.search(new ConstantScoreQuery(new MatchAllDocsQuery()), collector);
            return collector.getBitSet().toLongArray();
        }


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


    public synchronized ScoreHit[] simsearch(int fp[], int top, float cutoff) throws IOException
    {
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


    protected void close() throws IOException
    {
        try
        {
            if(searcher != null)
                searcher.getIndexReader().close();

            if(folder != null)
                folder.close();
        }
        finally
        {
            searcher = null;
            folder = null;
            idTable = null;
            sizeTable = null;
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
