package cz.iocb.orchem.load;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;



public class GeneralCompoundLoader
{
    public static void main(String[] args) throws Exception
    {
        if(args.length != 2)
        {
            System.err.println("wrong number of parameters");
            System.exit(1);
        }


        ConfigurationProperties properties = new ConfigurationProperties(args[0]);

        String pgHost = properties.getProperty("postgres.host");
        int pgPort = properties.getIntProperty("postgres.port");
        String pgUserName = properties.getProperty("postgres.username");
        String pgPassword = properties.getProperty("postgres.password");
        String pgDatabase = properties.getProperty("postgres.database");

        String idTag = properties.getProperty("sdf.idtag");
        String idPrefix = properties.getProperty("sdf.idprefix", "");

        String pgUrl = "jdbc:postgresql://" + pgHost + ":" + pgPort + "/" + pgDatabase;


        try(Connection connection = DriverManager.getConnection(pgUrl, pgUserName, pgPassword))
        {
            connection.setAutoCommit(false);

            try
            {
                CompoundLoader loader = new CompoundLoader(connection, idTag, idPrefix);
                loader.loadDirectory(new File(args[1]));

                try(Statement statement = connection.createStatement())
                {
                    statement.execute("select \"orchem_sync_data\"()");
                }

                connection.commit();
            }
            catch(Throwable e)
            {
                connection.rollback();
                throw e;
            }
        }
    }
}
