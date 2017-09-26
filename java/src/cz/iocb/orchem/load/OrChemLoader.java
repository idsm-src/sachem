package cz.iocb.orchem.load;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.BitSet;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.orchem.fingerprint.OrchemExtendedFingerprinter;
import cz.iocb.orchem.fingerprint.OrchemFingerprinter;
import cz.iocb.orchem.search.OrchemMoleculeBuilder;
import cz.iocb.orchem.shared.MoleculeCounts;
import cz.iocb.orchem.shared.MoleculeCreator;



public class OrChemLoader
{
    // bitmap index parameters
    private static final int fpSize = new OrchemFingerprinter().getSize();
    private static final int fpOffset = 1;

    // fingerprinter
    private static final IFingerprinter fingerprinter = new OrchemExtendedFingerprinter();

    // table names
    private static final String compoundsTable = "compounds";
    private static final String fingerprintTable = "fingerprint_orchem";
    private static final String countsTable = "molecule_counts";
    private static final String moleculesTable = "molecules";
    private static final String indexTable = "fingerprint_orchem_index";

    private static final int batchSize = 10000;


    public static void main(String[] args) throws Exception
    {
        int insertCount = 0;

        try (Connection insertConnection = ConnectionPool.getConnection())
        {
            try (PreparedStatement fingerprintInsert = insertConnection
                    .prepareStatement("insert into " + fingerprintTable + " (id, bit_count, fp) values (?,?,?)"))
            {
                try (PreparedStatement countsInsert = insertConnection.prepareStatement("insert into " + countsTable
                        + " (seqid, molTripleBondCount, molSCount, molOCount, molNCount, molFCount, molClCount, molBrCount, molICount, molCCount, molPCount) values (?,?,?,?,?,?,?,?,?,?,?)"))
                {
                    try (PreparedStatement moleculeInsert = insertConnection.prepareStatement(
                            "insert into " + moleculesTable + " (seqid, id, atoms, bonds) values (?,?,?,?)"))
                    {

                        try (Connection selectConnection = ConnectionPool.getConnection())
                        {
                            selectConnection.setAutoCommit(false);

                            // create empty bitmap index
                            int moleculeCount = 0;

                            try (Statement statement = selectConnection.createStatement(ResultSet.TYPE_FORWARD_ONLY,
                                    ResultSet.CONCUR_READ_ONLY))
                            {
                                try (ResultSet result = statement
                                        .executeQuery("select count(*) from " + compoundsTable))
                                {
                                    result.next();
                                    moleculeCount = result.getInt(1);
                                }
                            }

                            BitSet[] bitmasks = new BitSet[fpSize - fpOffset];

                            for(int i = 0; i < fpSize - fpOffset; i++)
                                bitmasks[i] = new BitSet(moleculeCount);


                            // load compounds
                            try (PreparedStatement selectStatement = selectConnection.prepareStatement(
                                    "SELECT id, molfile FROM " + compoundsTable + " ORDER BY id",
                                    ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY))
                            {
                                selectStatement.setFetchSize(batchSize);

                                try (ResultSet rs = selectStatement.executeQuery())
                                {
                                    while(rs.next())
                                    {
                                        int id = rs.getInt(1);
                                        String molfile = rs.getString(2);
                                        IAtomContainer readMolecule = MoleculeCreator.getMoleculeFromMolfile(molfile);
                                        MoleculeCreator.configureMolecule(readMolecule);

                                        // insert fingerprint
                                        BitSet fp = fingerprinter.getBitFingerprint(readMolecule).asBitSet();
                                        long[] words = fp.toLongArray();
                                        Object[] array = new Object[words.length];

                                        for(int i = 0; i < words.length; i++)
                                            array[i] = new Long(words[i]);

                                        fingerprintInsert.setInt(1, id);
                                        fingerprintInsert.setInt(2, fp.cardinality());
                                        fingerprintInsert.setArray(3, insertConnection.createArrayOf("bigint", array));
                                        fingerprintInsert.addBatch();

                                        // set bitmap indexes
                                        for(int idx = fp.nextSetBit(fpOffset); idx >= 0
                                                && idx < fpSize; idx = fp.nextSetBit(idx + 1))
                                            bitmasks[idx - fpOffset].set(insertCount);

                                        // insert molecule couts
                                        MoleculeCounts counts = new MoleculeCounts(readMolecule);
                                        countsInsert.setInt(1, insertCount);
                                        countsInsert.setShort(2, counts.molTripleBondCount);
                                        countsInsert.setShort(3, counts.molSCount);
                                        countsInsert.setShort(4, counts.molOCount);
                                        countsInsert.setShort(5, counts.molNCount);
                                        countsInsert.setShort(6, counts.molFCount);
                                        countsInsert.setShort(7, counts.molClCount);
                                        countsInsert.setShort(8, counts.molBrCount);
                                        countsInsert.setShort(9, counts.molICount);
                                        countsInsert.setShort(10, counts.molCCount);
                                        countsInsert.setShort(11, counts.molPCount);
                                        countsInsert.addBatch();

                                        // insert molecule binary representation
                                        OrchemMoleculeBuilder builder = new OrchemMoleculeBuilder(readMolecule);
                                        moleculeInsert.setInt(1, insertCount);
                                        moleculeInsert.setInt(2, id);
                                        moleculeInsert.setBytes(3, builder.atomsAsBytes());
                                        moleculeInsert.setBytes(4, builder.bondsAsBytes());
                                        moleculeInsert.addBatch();


                                        if(++insertCount % batchSize == 0)
                                        {
                                            fingerprintInsert.executeBatch();
                                            countsInsert.executeBatch();
                                            moleculeInsert.executeBatch();
                                        }
                                    }

                                    if(insertCount % batchSize != 0)
                                    {
                                        fingerprintInsert.executeBatch();
                                        countsInsert.executeBatch();
                                        moleculeInsert.executeBatch();
                                    }
                                }
                            }


                            // insert bitmap indexes
                            try (PreparedStatement insertStatement = insertConnection
                                    .prepareStatement("insert into " + indexTable + " (idx, bitmap) values (?,?)"))
                            {
                                for(int idx = 0; idx < bitmasks.length; idx++)
                                {
                                    insertStatement.setShort(1, (short) idx);
                                    insertStatement.setBytes(2, bitmasks[idx].toByteArray());
                                    insertStatement.addBatch();
                                }

                                insertStatement.executeBatch();
                            }
                        }
                    }
                }
            }
        }
    }
}
