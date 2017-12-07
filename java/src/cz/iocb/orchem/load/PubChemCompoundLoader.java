package cz.iocb.orchem.load;

import java.io.File;
import java.io.FileInputStream;
import java.sql.Connection;
import java.util.Properties;



public class PubChemCompoundLoader extends CompoundLoader
{
    public PubChemCompoundLoader(Connection connection)
    {
        super(connection, "> <PUBCHEM_COMPOUND_CID>", "");
    }


    public static void main(String[] args) throws Exception
    {
        if(args.length != 1)
        {
            System.err.println("missing directory name argument");
            System.exit(1);
        }


        String filename = System.getProperty("user.home") + "/.sachem/pubchem-datasource.properties";
        Properties properties = new Properties();

        try (FileInputStream stream = new FileInputStream(filename))
        {
            properties.load(stream);
        }

        ConnectionPool connectionPool = new ConnectionPool(properties);
        Connection connection = connectionPool.getConnection();
        connection.setAutoCommit(false);

        try
        {
            PubChemCompoundLoader loader = new PubChemCompoundLoader(connection);
            loader.loadDirectory(new File(args[0]));
            connection.commit();
        }
        catch (Throwable e)
        {
            connection.rollback();
        }
    }
}
