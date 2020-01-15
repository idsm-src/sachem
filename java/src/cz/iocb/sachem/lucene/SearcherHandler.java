package cz.iocb.sachem.lucene;

import static cz.iocb.sachem.lucene.Settings.idFieldName;
import java.io.IOException;
import java.util.Collections;
import java.util.Set;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;



public class SearcherHandler
{
    private static final Set<String> idSet = Collections.singleton(idFieldName);

    private IndexSearcher searcher;
    private Query query;
    private ScoreDoc after = null;

    public int[] ids;
    public float[] scores;


    public SearcherHandler(IndexSearcher searcher, Query query, int bufferSize)
    {
        this.searcher = searcher;
        this.query = query;
        this.ids = new int[bufferSize];
        this.scores = new float[bufferSize];
    }


    public int load(int limit) throws IOException
    {
        TopDocs hits = searcher.searchAfter(after, query, limit);
        ScoreDoc[] results = hits.scoreDocs;

        for(int i = 0; i < results.length; i++)
        {
            ids[i] = searcher.doc(results[i].doc, idSet).getField(idFieldName).numericValue().intValue();
            scores[i] = results[i].score;
        }

        if(results.length > 0)
            after = results[results.length - 1];

        return results.length;
    }
}
