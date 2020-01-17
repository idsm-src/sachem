package cz.iocb.sachem.load;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;



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
        String index = properties.getProperty("sachem.index");

        String idTag = properties.getProperty("sdf.idtag");
        String idPrefix = properties.getProperty("sdf.idprefix", "");

        String pgUrl = "jdbc:postgresql://" + pgHost + ":" + pgPort + "/" + pgDatabase;


        try(Connection connection = DriverManager.getConnection(pgUrl, pgUserName, pgPassword))
        {
            connection.setAutoCommit(false);

            try
            {
                CompoundLoader loader = new CompoundLoader(connection, index, idTag, idPrefix);
                loader.loadDirectory(new File(args[1]));

                try(PreparedStatement statement = connection.prepareStatement("select sachem.sync_data(?)"))
                {
                    statement.setString(1, index);
                    statement.execute();
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
