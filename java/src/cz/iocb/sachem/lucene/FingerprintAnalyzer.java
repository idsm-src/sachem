package cz.iocb.sachem.lucene;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.TokenStream;
import org.apache.lucene.analysis.Tokenizer;



public class FingerprintAnalyzer extends Analyzer
{
    @Override
    protected TokenStreamComponents createComponents(String fieldName)
    {
        Tokenizer source = new FingerprintTokenizer();
        return new TokenStreamComponents(source);
    }


    @Override
    protected TokenStream normalize(String fieldName, TokenStream in)
    {
        throw new UnsupportedOperationException();
    }
};
