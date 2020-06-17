package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.apache.lucene.index.DocValues;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.index.NumericDocValues;
import org.apache.lucene.search.Collector;
import org.apache.lucene.search.CollectorManager;
import org.apache.lucene.search.LeafCollector;
import org.apache.lucene.search.Scorable;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.ScoreMode;



public class SortedResultCollectorManager implements CollectorManager<Collector, SearchResult>
{
    private static final Comparator<ScoreDoc> comparator = new Comparator<ScoreDoc>()
    {
        @Override
        public int compare(ScoreDoc x, ScoreDoc y)
        {
            if(x.score == y.score)
                return x.doc - y.doc;

            return Float.compare(y.score, x.score);
        }
    };


    private final String name;
    private final List<ScoreDoc> results;


    SortedResultCollectorManager(String name)
    {
        this.name = name;
        results = Collections.synchronizedList(new ArrayList<ScoreDoc>());
    }


    @Override
    public Collector newCollector() throws IOException
    {
        return new Collector()
        {
            @Override
            public LeafCollector getLeafCollector(LeafReaderContext context) throws IOException
            {
                return new LeafCollector()
                {
                    NumericDocValues idField = DocValues.getNumeric(context.reader(), Settings.idFieldName);
                    Scorable scorer = null;

                    @Override
                    public void setScorer(Scorable scorer) throws IOException
                    {
                        this.scorer = scorer;
                    }

                    @Override
                    public void collect(int doc) throws IOException
                    {
                        idField.advanceExact(doc);

                        results.add(new ScoreDoc((int) idField.longValue(), scorer.score()));
                    }
                };
            }

            @Override
            public ScoreMode scoreMode()
            {
                return ScoreMode.COMPLETE;
            }
        };
    }


    @Override
    public SearchResult reduce(Collection<Collector> collectors) throws IOException
    {
        int ids[] = new int[results.size()];
        float scores[] = new float[results.size()];

        Collections.sort(results, comparator);

        for(int i = 0; i < results.size(); i++)
        {
            ScoreDoc node = results.get(i);

            ids[i] = node.doc;
            scores[i] = node.score;
        }

        return new SearchResult(name, results.size(), ids, scores);
    }
}
