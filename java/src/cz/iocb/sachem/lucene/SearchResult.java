package cz.iocb.sachem.lucene;



public class SearchResult
{
    public final int length;
    public final int[] ids;
    public final float[] scores;


    public SearchResult(int length, int[] ids, float[] scores)
    {
        this.length = length;
        this.ids = ids;
        this.scores = scores;
    }
}
