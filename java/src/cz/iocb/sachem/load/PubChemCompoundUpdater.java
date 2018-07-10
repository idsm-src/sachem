package cz.iocb.sachem.load;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
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



public class PubChemCompoundUpdater
{
    private static final int DAY = 24 * 60 * 60 * 1000;
    private static final String DAILY = "/Daily";
    private static final String WEEKLY = "/Weekly";


    public static void main(String[] args) throws Exception
    {
        if(args.length != 1)
        {
            System.err.println("wrong number of parameters");
            System.exit(1);
        }


        Date checkdate = new Date();
        ConfigurationProperties properties = new ConfigurationProperties(args[0]);

        String pgHost = properties.getProperty("postgres.host");
        int pgPort = properties.getIntProperty("postgres.port");
        String pgUserName = properties.getProperty("postgres.username");
        String pgPassword = properties.getProperty("postgres.password");
        String pgDatabase = properties.getProperty("postgres.database");
        boolean autoclean = properties.getBooleanProperty("sachem.autoclean");

        String ftpServer = properties.getProperty("ftp.server");
        int ftpPort = properties.getIntProperty("ftp.port");
        String ftpUserName = properties.getProperty("ftp.username");
        String ftpPassword = properties.getProperty("ftp.password");
        String ftpPath = properties.getProperty("ftp.path");

        String filePattern = properties.getProperty("sdf.pattern");
        String workdir = properties.getProperty("sdf.directory");
        String idTag = properties.getProperty("sdf.idtag");
        String idPrefix = properties.getProperty("sdf.idprefix");

        String baseDirectory = properties.getProperty("base.directory");
        String baseVersion = properties.getProperty("base.version");


        String pgUrl = "jdbc:postgresql://" + pgHost + ":" + pgPort + "/" + pgDatabase;

        try(Connection connection = DriverManager.getConnection(pgUrl, pgUserName, pgPassword))
        {
            if(autoclean)
            {
                try(Statement statement = connection.createStatement())
                {
                    statement.execute("select \"sachem_cleanup\"()");
                }
            }


            String loadedVersion = null;

            try(Statement statement = connection.createStatement())
            {
                try(ResultSet result = statement.executeQuery("select version from compound_stats where id = 0"))
                {
                    if(result.next())
                        loadedVersion = result.getString(1);
                }
            }


            LinkedList<String> updateList = new LinkedList<String>();
            String finalVersion = null;
            FTPClient ftpClient = new FTPClient();

            try
            {
                ftpClient.connect(ftpServer, ftpPort);
                ftpClient.login(ftpUserName, ftpPassword);
                ftpClient.setFileType(FTP.BINARY_FILE_TYPE);
                ftpClient.enterLocalPassiveMode();

                String lastVersion = loadedVersion != null ? loadedVersion : baseVersion;
                finalVersion = lastVersion;
                SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd");

                LinkedList<String> dailyList = new LinkedList<String>();
                LinkedList<String> weeklyList = new LinkedList<String>();


                FTPFile[] dailyUpdates = ftpClient.listFiles(ftpPath + "/Daily");

                if(dailyUpdates.length == 0)
                    throw new Exception("ftp directory is empty");

                String oldestDailyUpdate = null;
                boolean coveredByDaily = false;

                for(FTPFile update : dailyUpdates)
                {
                    String name = update.getName();

                    if(name.compareTo(lastVersion) == 0)
                    {
                        coveredByDaily = true;
                    }
                    if(name.compareTo(lastVersion) > 0)
                    {
                        dailyList.add(DAILY + "/" + name);

                        if(oldestDailyUpdate == null || name.compareTo(oldestDailyUpdate) < 0)
                            oldestDailyUpdate = name;

                        if(name.compareTo(finalVersion) > 0)
                            finalVersion = name;
                    }
                }


                if(!coveredByDaily)
                {
                    if(oldestDailyUpdate == null)
                        throw new Exception("inconsistent server daily data");

                    FTPFile[] weeklyUpdates = ftpClient.listFiles(ftpPath + "/Weekly");

                    if(weeklyUpdates.length == 0)
                        throw new Exception("ftp directory is empty");

                    String limit = format.format(new Date(format.parse(oldestDailyUpdate).getTime() + 6 * DAY));
                    String oldestWeeklyUpdate = null;

                    for(FTPFile update : weeklyUpdates)
                    {
                        String name = update.getName();

                        if(name.compareTo(lastVersion) > 0 && name.compareTo(limit) < 0)
                        {
                            weeklyList.add(WEEKLY + "/" + name);

                            if(oldestWeeklyUpdate == null || name.compareTo(oldestWeeklyUpdate) < 0)
                                oldestWeeklyUpdate = name;
                        }
                    }

                    if(oldestWeeklyUpdate == null)
                        throw new Exception("inconsistent server weekly data");

                    String firstUncovered = format
                            .format(new Date(format.parse(oldestWeeklyUpdate).getTime() - 7 * DAY));

                    if(firstUncovered.compareTo(lastVersion) > 0)
                        throw new Exception("database version is too old");
                }

                updateList.addAll(weeklyList);
                updateList.addAll(dailyList);


                for(String update : updateList)
                {
                    System.out.println("download " + update);

                    String basePath = workdir + update;
                    String sdfPath = basePath + "/" + "SDF";
                    File directory = new File(sdfPath);
                    directory.mkdirs();

                    try(FileOutputStream output = new FileOutputStream(basePath + "/killed-CIDs"))
                    {
                        ftpClient.retrieveFile(ftpPath + update + "/killed-CIDs", output);
                    }


                    FTPFile[] sdfUpdates = ftpClient.listFiles(ftpPath + update + "/SDF");

                    for(FTPFile sdfUpdate : sdfUpdates)
                    {
                        String name = sdfUpdate.getName();

                        if(name.matches(filePattern))
                        {
                            try(FileOutputStream output = new FileOutputStream(sdfPath + "/" + name))
                            {
                                System.out.println("  " + name);

                                if(!ftpClient.retrieveFile(ftpPath + update + "/SDF/" + name, output))
                                    throw new Exception("cannot download: " + ftpPath + update + "/SDF/" + name);
                            }
                        }
                    }
                }
            }
            finally
            {
                ftpClient.logout();
                ftpClient.disconnect();
            }


            if(loadedVersion != null && updateList.isEmpty())
            {
                try(PreparedStatement statement = connection
                        .prepareStatement("update compound_stats set checkdate=? where id = 0"))
                {
                    statement.setTimestamp(1, new Timestamp(checkdate.getTime()));
                    statement.executeUpdate();
                }

                return;
            }


            connection.setAutoCommit(false);

            try
            {
                CompoundLoader loader = new CompoundLoader(connection, idTag, idPrefix);

                if(loadedVersion == null)
                {
                    System.out.println("load " + baseDirectory);
                    loader.loadDirectory(new File(baseDirectory));
                }

                for(String update : updateList)
                {
                    try(PreparedStatement deleteStatement = connection
                            .prepareStatement("delete from compounds where id = ?"))
                    {
                        System.out.println("load " + workdir + update);

                        try(BufferedReader removeList = new BufferedReader(
                                new FileReader(workdir + update + "/killed-CIDs")))
                        {
                            int count = 0;
                            String line;

                            while((line = removeList.readLine()) != null)
                            {
                                count++;
                                deleteStatement.setInt(1, Integer.parseInt(line));
                                deleteStatement.addBatch();

                                if(count % 10000 == 0)
                                    deleteStatement.executeBatch();
                            }

                            if(count % 10000 != 0)
                                deleteStatement.executeBatch();
                        }
                    }

                    loader.loadDirectory(new File(workdir + update + "/SDF"), false);
                }

                try(Statement statement = connection.createStatement())
                {
                    statement.execute("select \"sachem_sync_data\"()");
                }

                try(PreparedStatement statement = connection
                        .prepareStatement("insert into compound_stats (id,version,size,checkdate) "
                                + "values (0,?,(select count(*) from compounds),?) on conflict (id) do update set "
                                + "version=EXCLUDED.version, size=EXCLUDED.size, checkdate=EXCLUDED.checkdate"))
                {
                    statement.setString(1, finalVersion);
                    statement.setTimestamp(2, new Timestamp(checkdate.getTime()));
                    statement.executeUpdate();
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
