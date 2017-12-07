package cz.iocb.orchem.load;

import java.io.File;
import java.io.FileOutputStream;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Statement;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.LinkedList;
import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;



public class CompoundUpdater
{
    public static void main(String[] args) throws Exception
    {
        if(args.length != 1)
        {
            System.err.println("usage: sachem-updatedb configfile.properties");
            System.exit(1);
        }


        ConfigurationProperties properties = new ConfigurationProperties(args[0]);

        String pgHost = properties.getProperty("postgres.host");
        int pgPort = properties.getIntProperty("postgres.port");
        String pgUserName = properties.getProperty("postgres.username");
        String pgPassword = properties.getProperty("postgres.password");
        String pgDatabase = properties.getProperty("postgres.database");

        String ftpServer = properties.getProperty("ftp.server");
        int ftpPort = properties.getIntProperty("ftp.port");
        String ftpUserName = properties.getProperty("ftp.username");
        String ftpPassword = properties.getProperty("ftp.password");
        String ftpPath = properties.getProperty("ftp.path");

        String filePattern = properties.getProperty("sdf.pattern");
        String workdir = properties.getProperty("sdf.directory");
        String idTag = properties.getProperty("sdf.idtag");
        String idPrefix = properties.getProperty("sdf.idprefix", "");


        String pgUrl = "jdbc:postgresql://" + pgHost + ":" + pgPort + "/" + pgDatabase;
        String path = workdir + "/" + new SimpleDateFormat("yyyy-MM-dd_HH:mm:ss").format(new Date());
        File directory = new File(path);

        LinkedList<FTPFile> sdfFiles = new LinkedList<FTPFile>();
        boolean hasNewItem = false;


        try (Connection connection = DriverManager.getConnection(pgUrl, pgUserName, pgPassword))
        {
            FTPClient ftpClient = new FTPClient();

            try
            {
                try (PreparedStatement statement = connection
                        .prepareStatement("select id from compound_sources where name=? and size=? and timestamp=?"))
                {
                    ftpClient.connect(ftpServer, ftpPort);
                    ftpClient.login(ftpUserName, ftpPassword);
                    ftpClient.setFileType(FTP.BINARY_FILE_TYPE);

                    for(FTPFile file : ftpClient.listFiles(ftpPath))
                    {
                        if(!file.isDirectory() && file.getName().matches(filePattern))
                        {
                            statement.setString(1, file.getName());
                            statement.setLong(2, file.getSize());
                            statement.setTimestamp(3, new Timestamp(file.getTimestamp().getTimeInMillis()));

                            sdfFiles.add(file);

                            try (ResultSet result = statement.executeQuery())
                            {
                                if(!result.next())
                                    hasNewItem = true;
                            }
                        }
                    }
                }

                if(hasNewItem)
                {
                    directory.mkdirs();

                    for(FTPFile sdfFile : sdfFiles)
                    {
                        String name = sdfFile.getName();

                        try (FileOutputStream output = new FileOutputStream(path + "/" + name))
                        {
                            if(!ftpClient.retrieveFile(ftpPath + "/" + name, output))
                                throw new Exception("cannot download: " + ftpPath + "/" + name);
                        }
                    }
                }
            }
            finally
            {
                ftpClient.logout();
                ftpClient.disconnect();
            }


            if(!hasNewItem)
                return;


            connection.setAutoCommit(false);

            try
            {
                CompoundLoader loader = new CompoundLoader(connection, idTag, idPrefix);
                loader.loadDirectory(directory);

                try (Statement statement = connection.createStatement())
                {
                    statement.execute("delete from compound_sources");
                }

                try (PreparedStatement statement = connection
                        .prepareStatement("insert into compound_sources (name, size, timestamp) values(?,?,?)"))
                {
                    for(FTPFile sdfFile : sdfFiles)
                    {
                        statement.setString(1, sdfFile.getName());
                        statement.setLong(2, sdfFile.getSize());
                        statement.setTimestamp(3, new Timestamp(sdfFile.getTimestamp().getTimeInMillis()));
                        statement.addBatch();
                    }

                    statement.executeBatch();
                }

                try (Statement statement = connection.createStatement())
                {
                    statement.execute("select \"orchem_sync_data\"()");
                }

                connection.commit();
            }
            catch (Throwable e)
            {
                connection.rollback();
                throw e;
            }
        }
    }
}
