package cz.iocb.orchem.search;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import cz.iocb.orchem.fingerprint.OrchemExtendedFingerprinter;
import cz.iocb.orchem.shared.MoleculeCreator;



public class OrchemSimilaritySearch
{
    private static final ThreadLocal<OrchemExtendedFingerprinter> fingerPrinter = new ThreadLocal<OrchemExtendedFingerprinter>()
    {
        @Override
        protected OrchemExtendedFingerprinter initialValue()
        {
            return new OrchemExtendedFingerprinter();
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
            molecule = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(query);
        else
            throw new CDKException("unsupported format");

        MoleculeCreator.configureMolecule(molecule);

        BitSet fp = fingerPrinter.get().getFingerprint(molecule, 10000);
        return fp.toLongArray();
    }
}
