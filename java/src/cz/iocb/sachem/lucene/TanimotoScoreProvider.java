package cz.iocb.sachem.lucene;

import java.io.IOException;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.queries.CustomScoreProvider;



public class TanimotoScoreProvider extends CustomScoreProvider
{
    private final Lucene lucene;
    private final int querySize;
    private final int docBase;


    public TanimotoScoreProvider(LeafReaderContext context, Lucene lucene, int querySize)
    {
        super(context);
        this.lucene = lucene;
        this.querySize = querySize;
        this.docBase = context.docBase;
    }


    @Override
    public float customScore(int doc, float sharedSize, float valSrcScores[]) throws IOException
    {
        int targetSize = lucene.getMoleculeFpSize(docBase + doc);

        return sharedSize / (querySize + targetSize - sharedSize);
    }
}
