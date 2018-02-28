package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.io.Reader;



public class FingerprintReader extends Reader
{
    int[] fingerprint;


    public FingerprintReader(int[] fp)
    {
        fingerprint = fp;
    }


    @Override
    public void close() throws IOException
    {
    }


    @Override
    public int read(char[] cbuf, int off, int len) throws IOException
    {
        throw new UnsupportedOperationException();
    }
}
