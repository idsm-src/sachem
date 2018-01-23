package cz.iocb.orchem.search;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.io.FormatFactory;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.formats.MDLV2000Format;
import org.openscience.cdk.io.formats.MDLV3000Format;
import org.openscience.cdk.io.formats.RGroupQueryFormat;
import org.openscience.cdk.io.formats.SMILESFormat;



public enum QueryFormat
{
    UNSPECIFIED, SMILES, MOLFILE, RGROUP;


    public static QueryFormat detect(String data) throws IOException, CDKException
    {
        IChemFormat format = new FormatFactory().guessFormat(new ByteArrayInputStream(data.getBytes()));

        if(format instanceof RGroupQueryFormat)
            return RGROUP;
        else if(format instanceof MDLV2000Format || format instanceof MDLV3000Format)
            return MOLFILE;
        else if(format instanceof SMILESFormat)
            return SMILES;
        else if(format == null && (data.contains("\n") || data.contains("\r")))
            return MOLFILE;
        else if(format == null && !data.contains("\n") && !data.contains("\r"))
            return SMILES;
        else
            throw new CDKException("unsupported format");
    }
}
