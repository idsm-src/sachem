package cz.iocb.sachem.lucene;

import java.io.IOException;
import org.apache.lucene.analysis.Tokenizer;
import org.apache.lucene.analysis.tokenattributes.CharTermAttribute;
import org.apache.lucene.analysis.tokenattributes.PositionIncrementAttribute;



public class FingerprintTokenizer extends Tokenizer
{
    private static final char[] b64str = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
            'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
            'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4',
            '5', '6', '7', '8', '9', '+', '/' };

    private final CharTermAttribute charTermAttribute = addAttribute(CharTermAttribute.class);
    private final PositionIncrementAttribute posIncrAtt = addAttribute(PositionIncrementAttribute.class);
    private final char[] buffer = new char[6];
    private int[] fingerprint;
    private int position = 0;


    public String bitAsString(int bit)
    {
        for(int i = 0; i < 6; i++)
            buffer[i] = b64str[bit >>> 6 * i & 0x3f];

        return new String(buffer);
    }


    @Override
    public boolean incrementToken() throws IOException
    {
        charTermAttribute.setEmpty();

        if(position == fingerprint.length)
            return false;

        int bit = fingerprint[position++];

        charTermAttribute.append(bitAsString(bit));
        posIncrAtt.setPositionIncrement(1);
        return true;
    }


    @Override
    public void reset() throws IOException
    {
        super.reset();

        fingerprint = ((FingerprintReader) input).fingerprint;
        position = 0;
    }
}
