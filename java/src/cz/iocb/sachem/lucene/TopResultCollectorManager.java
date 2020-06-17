package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import org.apache.lucene.index.DocValues;
import org.apache.lucene.index.IndexReaderContext;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.index.NumericDocValues;
import org.apache.lucene.search.Collector;
import org.apache.lucene.search.CollectorManager;
import org.apache.lucene.search.LeafCollector;
import org.apache.lucene.search.Scorable;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.ScoreMode;



public class TopResultCollectorManager implements CollectorManager<Collector, SearchResult>
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
    private int limit;
    private int hits = 0;
    private ScoreDoc[] heap;
    private ArrayList<Integer> timeouted = new ArrayList<Integer>();


    TopResultCollectorManager(String name, int limit)
    {
        this.name = name;
        this.limit = limit;
    }


    @Override
    public Collector newCollector() throws IOException
    {
        return new Collector()
        {
            @Override
            public LeafCollector getLeafCollector(LeafReaderContext context) throws IOException
            {
                if(heap == null)
                {
                    synchronized(TopResultCollectorManager.this)
                    {
                        if(heap == null)
                        {
                            IndexReaderContext parent = context;

                            while(parent.parent != null)
                                parent = parent.parent;

                            limit = Math.min(limit, parent.reader().maxDoc());
                            heap = new ScoreDoc[limit];

                            System.out.println("limit = " + limit);
                        }
                    }
                }


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
                        float score = scorer.score();

                        idField.advanceExact(doc);
                        int id = (int) idField.longValue();

                        synchronized(TopResultCollectorManager.this)
                        {
                            if(score == Float.NEGATIVE_INFINITY)
                            {
                                timeouted.add(id);
                            }
                            if(hits < limit)
                            {
                                heap[hits] = new ScoreDoc(id, score);

                                int j = hits++;

                                while(j != 0)
                                {
                                    int p = (j - 1) / 2;

                                    if(comparator.compare(heap[j], heap[p]) < 0)
                                        break;

                                    ScoreDoc node = heap[j];
                                    heap[j] = heap[p];
                                    heap[p] = node;

                                    j = p;
                                }
                            }
                            else if(score > heap[0].score || score == heap[0].score && id < heap[0].doc)
                            {
                                ScoreDoc node = heap[0];

                                node.doc = id;
                                node.score = score;

                                int i = 0;
                                int j = ((i + 1) << 1) - 1;
                                int k = j + 1;

                                if(k < limit && comparator.compare(heap[k], heap[j]) > 0)
                                    j = k;

                                while(j < limit && comparator.compare(heap[j], node) > 0)
                                {
                                    heap[i] = heap[j];

                                    i = j;
                                    j = ((i + 1) << 1) - 1;
                                    k = j + 1;

                                    if(k < limit && comparator.compare(heap[k], heap[j]) > 0)
                                        j = k;
                                }

                                heap[i] = node;
                            }
                        }
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
        int ids[] = new int[hits + timeouted.size()];
        float scores[] = new float[hits + timeouted.size()];

        Arrays.sort(heap, 0, hits, comparator);
        Collections.sort(timeouted);

        for(int i = 0; i < hits; i++)
        {
            ids[i] = heap[i].doc;
            scores[i] = heap[i].score;
        }

        for(int i = 0; i < timeouted.size(); i++)
        {
            ids[hits + i] = timeouted.get(i);
            scores[hits + i] = Float.NEGATIVE_INFINITY;
        }

        return new SearchResult(name, hits + timeouted.size(), ids, scores);
    }
}
