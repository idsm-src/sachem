package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.util.BitSet;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.search.ScoreMode;
import org.apache.lucene.search.SimpleCollector;



public class BitSetCollector extends SimpleCollector
{
    private final Lucene lucene;
    private final BitSet bitset;
    private int docBase;


    BitSetCollector(Lucene lucene, int maxMoleculeId)
    {
        this.lucene = lucene;
        this.bitset = new BitSet(maxMoleculeId);
    }


    @Override
    protected void doSetNextReader(LeafReaderContext context) throws IOException
    {
        docBase = context.docBase;
    }


    @Override
    public void collect(int docId) throws IOException
    {
        bitset.set(lucene.getMoleculeId(docBase + docId));
    }


    @Override
    public ScoreMode scoreMode()
    {
        return ScoreMode.COMPLETE_NO_SCORES;
    }


    public BitSet getBitSet()
    {
        return bitset;
    }
}
