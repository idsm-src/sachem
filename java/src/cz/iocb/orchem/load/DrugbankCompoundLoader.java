package cz.iocb.orchem.load;

import java.io.File;
import java.io.FileInputStream;
import java.sql.Connection;
import java.util.Properties;



public class DrugbankCompoundLoader extends CompoundLoader
{
    public DrugbankCompoundLoader(Connection connection)
    {
        super(connection, "> <DRUGBANK_ID>", "DB");
    }


    public static void main(String[] args) throws Exception
    {
        if(args.length != 1)
        {
            System.err.println("missing directory name argument");
            System.exit(1);
        }


        String filename = System.getProperty("user.home") + "/.sachem/drugbank-datasource.properties";
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
            DrugbankCompoundLoader loader = new DrugbankCompoundLoader(connection);
            loader.loadDirectory(new File(args[0]));
            connection.commit();
        }
        catch (Throwable e)
        {
            connection.rollback();
        }
    }
}
