package cz.iocb.orchem.load;

import java.io.File;
import java.io.FileInputStream;
import java.sql.Connection;
import java.util.Properties;



public class ChebiCompoundLoader extends CompoundLoader
{
    public ChebiCompoundLoader(Connection connection)
    {
        super(connection, "> <ChEBI ID>", "CHEBI:");
    }


    public static void main(String[] args) throws Exception
    {
        if(args.length != 1)
        {
            System.err.println("missing directory name argument");
            System.exit(1);
        }


        String filename = System.getProperty("user.home") + "/.sachem/chebi-datasource.properties";
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
            ChebiCompoundLoader loader = new ChebiCompoundLoader(connection);
            loader.loadDirectory(new File(args[0]));
            connection.commit();
        }
        catch (Throwable e)
        {
            connection.rollback();
        }
    }
}
