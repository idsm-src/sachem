package cz.iocb.sachem.search;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.concurrent.TimeoutException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.sachem.molecule.AromaticityMode;
import cz.iocb.sachem.molecule.BinaryMoleculeBuilder;
import cz.iocb.sachem.molecule.MoleculeCreator;
import cz.iocb.sachem.molecule.TautomerMode;



public class SimilaritySearch
{
    public static byte[] getQueryData(byte[] queryArray, int queryTypeIdx)
            throws IOException, CDKException, TimeoutException
    {
        String query = new String(queryArray, StandardCharsets.ISO_8859_1);

        List<IAtomContainer> queries = MoleculeCreator.translateQuery(query, AromaticityMode.PRESERVE,
                TautomerMode.IGNORE);

        BinaryMoleculeBuilder builder = new BinaryMoleculeBuilder(queries.get(0));
        return builder.asBytes(false);
    }
}
