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
import cz.iocb.orchem.isomorphism.IsomorphismSort;
import cz.iocb.orchem.search.OrchemMoleculeBuilder;
import cz.iocb.orchem.shared.MoleculeCounts;
import cz.iocb.orchem.shared.MoleculeCreator;



public class OrChemLoader
{
    static class Item
    {
        int id;
        int seqid;
        String molfile;

        MoleculeCounts counts;

        Object[] similarityFingerprint;
        Object[] substructureFingerprint;
        int similarityFingerprintCardinality;

        byte[] atoms;
        byte[] bonds;
    }


    // table names
    private static final String compoundsTable = "compounds";
    private static final String moleculesTable = "orchem_molecules";
    private static final String moleculeCountsTable = "orchem_molecule_counts";
    private static final String similarityFingerprintTable = "orchem_similarity_fingerprint";
    private static final String substructureFingerprintTable = "orchem_substructure_fingerprint";
    private static final String substructureFingerprintIndexTable = "orchem_substructure_fingerprint_index";

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
            try (PreparedStatement moleculesStatement = insertConnection
                    .prepareStatement("insert into " + moleculesTable + " (seqid, id, atoms, bonds) values (?,?,?,?)"))
            {
                try (PreparedStatement moleculeCountsStatement = insertConnection.prepareStatement("insert into "
                        + moleculeCountsTable
                        + " (seqid, molTripleBondCount, molSCount, molOCount, molNCount, molFCount, molClCount, molBrCount, molICount, molCCount, molPCount) values (?,?,?,?,?,?,?,?,?,?,?)"))
                {
                    try (PreparedStatement similarityFingerprintStatement = insertConnection.prepareStatement(
                            "insert into " + similarityFingerprintTable + " (id, bit_count, fp) values (?,?,?)"))
                    {
                        try (PreparedStatement substructureFingerprintStatement = insertConnection.prepareStatement(
                                "insert into " + substructureFingerprintTable + " (seqid, fp) values (?,?)"))
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
                                                // insert molecule binary representation
                                                moleculesStatement.setInt(1, item.seqid);
                                                moleculesStatement.setInt(2, item.id);
                                                moleculesStatement.setBytes(3, item.atoms);
                                                moleculesStatement.setBytes(4, item.bonds);
                                                moleculesStatement.addBatch();

                                                // insert molecule couts
                                                moleculeCountsStatement.setInt(1, item.seqid);
                                                moleculeCountsStatement.setShort(2, item.counts.molTripleBondCount);
                                                moleculeCountsStatement.setShort(3, item.counts.molSCount);
                                                moleculeCountsStatement.setShort(4, item.counts.molOCount);
                                                moleculeCountsStatement.setShort(5, item.counts.molNCount);
                                                moleculeCountsStatement.setShort(6, item.counts.molFCount);
                                                moleculeCountsStatement.setShort(7, item.counts.molClCount);
                                                moleculeCountsStatement.setShort(8, item.counts.molBrCount);
                                                moleculeCountsStatement.setShort(9, item.counts.molICount);
                                                moleculeCountsStatement.setShort(10, item.counts.molCCount);
                                                moleculeCountsStatement.setShort(11, item.counts.molPCount);
                                                moleculeCountsStatement.addBatch();

                                                // insert similarity fingerprint
                                                similarityFingerprintStatement.setInt(1, item.id);
                                                similarityFingerprintStatement.setInt(2,
                                                        item.similarityFingerprintCardinality);
                                                similarityFingerprintStatement.setArray(3, insertConnection
                                                        .createArrayOf("bigint", item.similarityFingerprint));
                                                similarityFingerprintStatement.addBatch();

                                                // insert substructure fingerprint
                                                substructureFingerprintStatement.setInt(1, item.seqid);
                                                substructureFingerprintStatement.setArray(2, insertConnection
                                                        .createArrayOf("smallint", item.substructureFingerprint));
                                                substructureFingerprintStatement.addBatch();
                                            }

                                            moleculesStatement.executeBatch();
                                            moleculeCountsStatement.executeBatch();
                                            similarityFingerprintStatement.executeBatch();
                                            substructureFingerprintStatement.executeBatch();
                                        }
                                    }
                                }


                                // insert bitmap indexes
                                try (PreparedStatement substructureFingerprintInsexStatement = insertConnection
                                        .prepareStatement("insert into " + substructureFingerprintIndexTable
                                                + " (idx, bitmap) values (?,?)"))
                                {
                                    for(int idx = 0; idx < bitmasks.length; idx++)
                                    {
                                        substructureFingerprintInsexStatement.setShort(1, (short) idx);
                                        substructureFingerprintInsexStatement.setBytes(2, bitmasks[idx].toByteArray());
                                        substructureFingerprintInsexStatement.addBatch();
                                    }

                                    substructureFingerprintInsexStatement.executeBatch();
                                }
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

                // calculate molecule binary representation
                readMolecule.setAtoms(IsomorphismSort.atomsByFrequency(readMolecule));
                OrchemMoleculeBuilder builder = new OrchemMoleculeBuilder(readMolecule);
                item.atoms = builder.atomsAsBytes();
                item.bonds = builder.bondsAsBytes();


                // calculate molecule couts
                item.counts = new MoleculeCounts(readMolecule);


                // calculate similarity fingerprint
                BitSet fp = fingerPrinter.get().getBitFingerprint(readMolecule).asBitSet();
                long[] words = fp.toLongArray();
                Object[] array = new Object[words.length];

                for(int i = 0; i < words.length; i++)
                    array[i] = new Long(words[i]);

                item.similarityFingerprint = array;
                item.similarityFingerprintCardinality = fp.cardinality();


                // calculate substructure fingerprint
                int size = 0;
                for(int i = fp.nextSetBit(fpOffset); i > 0 && i < fpSize; i = fp.nextSetBit(i + 1))
                    size++;

                item.substructureFingerprint = new Object[size];

                int position = 0;
                for(int i = fp.nextSetBit(fpOffset); i > 0 && i < fpSize; i = fp.nextSetBit(i + 1))
                    item.substructureFingerprint[position++] = i - fpOffset;


                // set bitmap indexes
                synchronized(OrChemLoader.class)
                {
                    for(int idx = fp.nextSetBit(fpOffset); idx >= 0 && idx < fpSize; idx = fp.nextSetBit(idx + 1))
                        bitmasks[idx - fpOffset].set(item.seqid);
                }
            }
            catch (Throwable e)
            {
                e.printStackTrace();
            }
        });
    }
}
