package cz.iocb.orchem.load;

import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.orchem.fingerprint.OrchemExtendedFingerprinter;
import cz.iocb.orchem.fingerprint.OrchemFingerprinter;
import cz.iocb.orchem.search.OrchemMoleculeBuilder;
import cz.iocb.orchem.shared.MoleculeCounts;
import cz.iocb.orchem.shared.MoleculeCreator;



public class OrChemLoader
{
    static class Item
    {
        int id;
        String molfile;

        int seqid;

        MoleculeCounts counts;

        Object[] fp;
        int cardinality;

        byte[] atoms;
        byte[] bonds;
    }


    // table names
    private static final String compoundsTable = "compounds";
    private static final String fingerprintTable = "fingerprint_orchem";
    private static final String countsTable = "molecule_counts";
    private static final String moleculesTable = "molecules";
    private static final String indexTable = "fingerprint_orchem_index";

    // bitmap index parameters
    private static final int fpSize = new OrchemFingerprinter().getSize();
    private static final int fpOffset = 1;

    // fingerprinter
    private static final ThreadLocal<OrchemExtendedFingerprinter> fingerPrinter = new ThreadLocal<OrchemExtendedFingerprinter>()
    {
        @Override
        protected OrchemExtendedFingerprinter initialValue()
        {
            return new OrchemExtendedFingerprinter();
        }
    };

    private static final int batchSize = 10000;


    public static void main(String[] args) throws Exception
    {
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
                                    int seqid = 0;

                                    while(true)
                                    {
                                        ArrayList<Item> items = new ArrayList<Item>(batchSize);

                                        while(items.size() < batchSize && rs.next())
                                        {
                                            Item item = new Item();

                                            item.id = rs.getInt(1);
                                            item.molfile = rs.getString(2);
                                            item.seqid = seqid++;

                                            items.add(item);
                                        }

                                        if(items.size() == 0)
                                            break;

                                        process(items, bitmasks);

                                        for(Item item : items)
                                        {
                                            // insert fingerprint
                                            fingerprintInsert.setInt(1, item.id);
                                            fingerprintInsert.setInt(2, item.cardinality);
                                            fingerprintInsert.setArray(3,
                                                    insertConnection.createArrayOf("bigint", item.fp));
                                            fingerprintInsert.addBatch();

                                            // insert molecule couts
                                            countsInsert.setInt(1, item.seqid);
                                            countsInsert.setShort(2, item.counts.molTripleBondCount);
                                            countsInsert.setShort(3, item.counts.molSCount);
                                            countsInsert.setShort(4, item.counts.molOCount);
                                            countsInsert.setShort(5, item.counts.molNCount);
                                            countsInsert.setShort(6, item.counts.molFCount);
                                            countsInsert.setShort(7, item.counts.molClCount);
                                            countsInsert.setShort(8, item.counts.molBrCount);
                                            countsInsert.setShort(9, item.counts.molICount);
                                            countsInsert.setShort(10, item.counts.molCCount);
                                            countsInsert.setShort(11, item.counts.molPCount);
                                            countsInsert.addBatch();

                                            // insert molecule binary representation
                                            moleculeInsert.setInt(1, item.seqid);
                                            moleculeInsert.setInt(2, item.id);
                                            moleculeInsert.setBytes(3, item.atoms);
                                            moleculeInsert.setBytes(4, item.bonds);
                                            moleculeInsert.addBatch();
                                        }

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


    private static void process(ArrayList<Item> items, BitSet[] bitmasks) throws CDKException, IOException
    {
        items.stream().parallel().forEach(item -> {
            try
            {
                IAtomContainer readMolecule = MoleculeCreator.getMoleculeFromMolfile(item.molfile);
                MoleculeCreator.configureMolecule(readMolecule);

                // insert fingerprint
                BitSet fp = fingerPrinter.get().getBitFingerprint(readMolecule).asBitSet();
                long[] words = fp.toLongArray();
                Object[] array = new Object[words.length];

                for(int i = 0; i < words.length; i++)
                    array[i] = new Long(words[i]);

                item.fp = array;
                item.cardinality = fp.cardinality();

                // set bitmap indexes
                for(int idx = fp.nextSetBit(fpOffset); idx >= 0 && idx < fpSize; idx = fp.nextSetBit(idx + 1))
                    bitmasks[idx - fpOffset].set(item.seqid);

                // insert molecule couts
                item.counts = new MoleculeCounts(readMolecule);

                // insert molecule binary representation
                OrchemMoleculeBuilder builder = new OrchemMoleculeBuilder(readMolecule);
                item.atoms = builder.atomsAsBytes();
                item.bonds = builder.bondsAsBytes();
            }
            catch (Throwable e)
            {
                e.printStackTrace();
            }
        });
    }
}
