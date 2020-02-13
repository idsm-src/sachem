package cz.iocb.sachem.lucene;

import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;



public class SearchResult
{
    public final int length;
    public final int[] ids;
    public final float[] scores;


    public SearchResult()
    {
        this.length = 0;
        this.ids = new int[0];
        this.scores = new float[0];
    }


    public SearchResult(int length, int[] ids, float[] scores)
    {
        this.length = length;
        this.ids = ids;
        this.scores = scores;
    }


    public SearchResult(TopDocs hits)
    {
        this.length = hits.scoreDocs.length;
        this.ids = new int[length];
        this.scores = new float[length];

        int i = 0;

        for(ScoreDoc hit : hits.scoreDocs)
        {
            ids[i] = hit.doc;
            scores[i] = hit.score;
            i++;
        }
    }
}
