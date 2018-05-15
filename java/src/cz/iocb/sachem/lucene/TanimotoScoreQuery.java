package cz.iocb.sachem.lucene;

import java.io.IOException;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.queries.CustomScoreProvider;
import org.apache.lucene.queries.CustomScoreQuery;
import org.apache.lucene.search.BooleanQuery;



public class TanimotoScoreQuery extends CustomScoreQuery
{
    private final Lucene lucene;
    private final int querySize;


    public TanimotoScoreQuery(Lucene lucene, BooleanQuery subQuery)
    {
        super(subQuery);
        this.lucene = lucene;
        querySize = subQuery.clauses().size();
    }


    @Override
    protected CustomScoreProvider getCustomScoreProvider(LeafReaderContext context) throws IOException
    {
        return new TanimotoScoreProvider(context, lucene, querySize);
    }
}
