package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import org.apache.lucene.index.DocValues;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.index.NumericDocValues;
import org.apache.lucene.search.Collector;
import org.apache.lucene.search.CollectorManager;
import org.apache.lucene.search.LeafCollector;
import org.apache.lucene.search.Scorable;
import org.apache.lucene.search.ScoreMode;
import cz.iocb.sachem.lucene.ResultCollectorManager.ResultCollector;



public class ResultCollectorManager implements CollectorManager<ResultCollector, SearchResult>
{
    private static final int bufferSize = 100000;


    private static class ResultBuffer
    {
        int[] ids = new int[bufferSize];
        float[] scores = new float[bufferSize];
        int possition = 0;
    }


    static class ResultCollector implements Collector
    {
        public final List<ResultBuffer> results = new LinkedList<ResultBuffer>();
        private ResultBuffer result = null;

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
                    if(result == null)
                    {
                        result = new ResultBuffer();
                        results.add(result);
                    }

                    idField.advanceExact(doc);

                    result.ids[result.possition] = (int) idField.longValue();
                    result.scores[result.possition] = scorer.score();
                    result.possition++;

                    if(result.possition == bufferSize)
                        result = null;
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


    @Override
    public ResultCollector newCollector() throws IOException
    {
        return new ResultCollector();
    }


    @Override
    public SearchResult reduce(Collection<ResultCollector> collectors) throws IOException
    {
        int length = 0;

        for(ResultCollector collector : collectors)
            for(ResultBuffer result : collector.results)
                length += result.possition;

        int[] ids = new int[length];
        float[] scores = new float[length];

        int position = 0;

        for(ResultCollector collector : collectors)
        {
            for(ResultBuffer result : collector.results)
            {
                for(int i = 0; i < result.possition; i++)
                {
                    ids[position] = result.ids[i];
                    scores[position] = result.scores[i];
                    position++;
                }
            }
        }

        return new SearchResult(length, ids, scores);
    }
}
