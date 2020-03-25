package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.util.Collections;
import java.util.Set;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;



public class SearchResult
{
    private static final Set<String> fieldsToLoad = Collections.singleton(Settings.idFieldName);

    public final String name;
    public final int length;
    public final int[] ids;
    public final float[] scores;


    public SearchResult()
    {
        this.name = "";
        this.length = 0;
        this.ids = new int[0];
        this.scores = new float[0];
    }


    public SearchResult(String name, int length, int[] ids, float[] scores)
    {
        this.name = name;
        this.length = length;
        this.ids = ids;
        this.scores = scores;
    }


    public SearchResult(String name, IndexSearcher searcher, TopDocs hits) throws IOException
    {
        this.name = name;
        this.length = hits.scoreDocs.length;
        this.ids = new int[length];
        this.scores = new float[length];

        int i = 0;

        for(ScoreDoc hit : hits.scoreDocs)
        {
            ids[i] = searcher.doc(hit.doc, fieldsToLoad).getField(Settings.idFieldName).numericValue().intValue();
            scores[i] = hit.score;
            i++;
        }
    }
}
