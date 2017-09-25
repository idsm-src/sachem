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



public class SimilaritySearch
{
    private static final String QUERY_TYPE_MOL = "MOL";
    private static final String QUERY_TYPE_SMILES = "SMILES";


    private static final ThreadLocal<OrchemExtendedFingerprinter> fingerPrinter = new ThreadLocal<OrchemExtendedFingerprinter>()
    {
        @Override
        protected OrchemExtendedFingerprinter initialValue()
        {
            return new OrchemExtendedFingerprinter();
        }
    };


    public static long[] getQueryData(byte[] queryArray, String queryType) throws CDKException, IOException
    {
        String query = new String(queryArray, StandardCharsets.ISO_8859_1);
        IAtomContainer molecule = null;

        if(queryType.equals(QUERY_TYPE_MOL))
            molecule = MoleculeCreator.getMoleculeFromMolfile(query);
        else if(queryType.equals(QUERY_TYPE_SMILES))
            molecule = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(query);
        else
            return null;

        MoleculeCreator.configureMolecule(molecule);

        BitSet fp = fingerPrinter.get().getBitFingerprint(molecule).asBitSet();
        return fp.toLongArray();
    }
}
