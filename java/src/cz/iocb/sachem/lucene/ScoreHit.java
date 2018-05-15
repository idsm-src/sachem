package cz.iocb.sachem.lucene;



public class ScoreHit
{
    public int id;
    public float score;


    public ScoreHit(int id, float score)
    {
        this.id = id;
        this.score = score;
    }
}
