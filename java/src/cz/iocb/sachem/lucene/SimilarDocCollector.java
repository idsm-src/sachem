package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.search.Scorer;
import org.apache.lucene.search.SimpleCollector;
import org.apache.lucene.util.PriorityQueue;



public class SimilarDocCollector extends SimpleCollector
{
    private final Lucene lucene;
    private final PriorityQueue<ScoreHit> queue;
    private final ArrayList<ScoreHit> list;
    private final int top;
    private final float cutoff;
    private Scorer scorer;
    private int docBase;


    SimilarDocCollector(Lucene lucene, int top, float cutoff)
    {
        this.lucene = lucene;
        this.top = top;
        this.cutoff = cutoff;

        if(top <= 0)
        {
            list = new ArrayList<ScoreHit>();
            queue = null;
        }
        else
        {
            list = null;
            queue = new PriorityQueue<ScoreHit>(top)
            {
                @Override
                protected boolean lessThan(ScoreHit a, ScoreHit b)
                {
                    return a.score < b.score;
                }
            };
        }
    }


    @Override
    protected void doSetNextReader(LeafReaderContext context) throws IOException
    {
        docBase = context.docBase;
    }


    @Override
    public void collect(int docId) throws IOException
    {
        float score = scorer.score();

        if(score < cutoff)
            return;


        int id = lucene.getMoleculeId(docBase + docId);

        if(list != null)
        {
            list.add(new ScoreHit(id, score));
        }
        else if(queue.size() < top)
        {
            queue.add(new ScoreHit(id, score));
        }
        else if(queue.top().score < score)
        {
            queue.top().id = id;
            queue.top().score = score;
            queue.updateTop();
        }
    }


    @Override
    public void setScorer(Scorer scorer) throws IOException
    {
        this.scorer = scorer;
    }


    @Override
    public boolean needsScores()
    {
        return true;
    }


    public ScoreHit[] getDocs()
    {
        if(list != null)
        {
            Collections.sort(list, new Comparator<ScoreHit>()
            {
                @Override
                public int compare(ScoreHit a, ScoreHit b)
                {
                    return Float.compare(b.score, a.score);
                }
            });

            return list.toArray(new ScoreHit[0]);
        }
        else
        {
            int size = queue.size();
            ScoreHit[] array = new ScoreHit[size];

            for(int i = size - 1; i >= 0; i--)
                array[i] = queue.pop();

            return array;
        }
    }
}
