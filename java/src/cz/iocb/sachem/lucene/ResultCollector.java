package cz.iocb.sachem.lucene;

import java.io.IOException;
import org.apache.lucene.index.DocValues;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.index.NumericDocValues;
import org.apache.lucene.search.Collector;
import org.apache.lucene.search.LeafCollector;
import org.apache.lucene.search.Scorable;
import org.apache.lucene.search.ScoreMode;



public class ResultCollector implements Collector
{
    private static int initLength = 100000;

    int[] ids = new int[initLength];
    float[] scores = new float[initLength];
    int possition = 0;


    @Override
    public LeafCollector getLeafCollector(LeafReaderContext context) throws IOException
    {
        return new LeafCollector()
        {
            NumericDocValues idField = DocValues.getNumeric(context.reader(), Settings.idFieldName);
            Scorable scorer = null;

            @Override
            public void collect(int doc) throws IOException
            {
                if(possition == ids.length)
                {
                    int newLength = 2 * ids.length;

                    int[] newIds = new int[newLength];
                    System.arraycopy(ids, 0, newIds, 0, ids.length);
                    ids = newIds;

                    float[] newScores = new float[newLength];
                    System.arraycopy(scores, 0, newScores, 0, scores.length);
                    scores = newScores;
                }

                idField.advanceExact(doc);

                ids[possition] = (int) idField.longValue();
                scores[possition] = scorer.score();

                possition++;
            }

            @Override
            public void setScorer(Scorable scorer) throws IOException
            {
                this.scorer = scorer;
            }
        };
    }


    @Override
    public ScoreMode scoreMode()
    {
        return ScoreMode.COMPLETE;
    }
}
