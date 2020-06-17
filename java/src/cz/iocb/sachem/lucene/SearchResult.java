package cz.iocb.sachem.lucene;



public class SearchResult
{
    public final String name;
    public final int length;
    public final int[] ids;
    public final float[] scores;


    public SearchResult(String name)
    {
        this.name = name;
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
}
