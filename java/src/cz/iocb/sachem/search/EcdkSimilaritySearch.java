package cz.iocb.sachem.search;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.sachem.shared.MoleculeCreator;



public class EcdkSimilaritySearch
{
    private static final ThreadLocal<ExtendedFingerprinter> fingerPrinter = new ThreadLocal<ExtendedFingerprinter>()
    {
        @Override
        protected ExtendedFingerprinter initialValue()
        {
            return new ExtendedFingerprinter();
        }
    };


    public static long[] getQueryData(byte[] queryArray, int queryTypeIdx) throws CDKException, IOException
    {
        QueryFormat queryType = QueryFormat.values()[queryTypeIdx];
        String query = new String(queryArray, StandardCharsets.ISO_8859_1);
        IAtomContainer molecule = null;

        if(queryType == QueryFormat.UNSPECIFIED)
            queryType = QueryFormat.detect(query);


        if(queryType == QueryFormat.MOLFILE)
            molecule = MoleculeCreator.getMoleculeFromMolfile(query);
        else if(queryType == QueryFormat.SMILES)
            molecule = MoleculeCreator.getMoleculeFromSmiles(query);
        else
            throw new CDKException("unsupported format");

        MoleculeCreator.configureMolecule(molecule);

        BitSet fp = fingerPrinter.get().getFingerprint(molecule);
        return fp.toLongArray();
    }
}
