package cz.iocb.sachem.search;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.sachem.isomorphism.IsomorphismSort;
import cz.iocb.sachem.shared.MoleculeCreator;



public class SimilaritySearch
{
    public static byte[] getQueryData(byte[] queryArray, int queryTypeIdx) throws IOException, CDKException
    {
        String query = new String(queryArray, StandardCharsets.ISO_8859_1);
        QueryFormat queryType = QueryFormat.values()[queryTypeIdx];

        if(queryType == QueryFormat.UNSPECIFIED)
            queryType = QueryFormat.detect(query);


        IAtomContainer molecule = null;

        if(queryType == QueryFormat.MOLFILE)
            molecule = MoleculeCreator.getMoleculeFromMolfile(query);
        else if(queryType == QueryFormat.SMILES)
            molecule = MoleculeCreator.getMoleculeFromSmiles(query);
        else
            throw new IllegalArgumentException("unsupported query format");


        MoleculeCreator.configureMolecule(molecule);

        molecule.setAtoms(IsomorphismSort.atomsByFrequency(molecule));
        SachemMoleculeBuilder builder = new SachemMoleculeBuilder(molecule);
        return builder.asBytes(false);
    }
}
