package cz.iocb.orchem.load;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.Arrays;
import java.util.Comparator;
import java.util.zip.GZIPInputStream;



public class CompoundLoader
{
    private static final int batchSize = 1000;
    private final Connection connection;
    private final String idTag;
    private final String idPrefix;


    protected CompoundLoader(Connection connection, String idTag, String idPrefix)
    {
        this.connection = connection;
        this.idTag = idTag;
        this.idPrefix = idPrefix;
    }


    private void parse(InputStream inputStream, int version) throws Exception
    {
        Reader decoder = new InputStreamReader(inputStream, "US-ASCII");
        BufferedReader reader = new BufferedReader(decoder);
        String line;

        try (PreparedStatement insertStatement = connection
                .prepareStatement("insert into compounds (id, version, molfile) values (?,?,?) "
                        + "on conflict (id) do update set version=EXCLUDED.version, molfile=EXCLUDED.molfile"))
        {
            int count = 0;

            while((line = reader.readLine()) != null)
            {
                Integer id = null;
                StringBuilder sdfBuilder = new StringBuilder();

                while(!line.startsWith(">") && !line.equals("$$$$") && line != null)
                {
                    sdfBuilder.append(line);
                    sdfBuilder.append('\n');
                    line = reader.readLine();
                }

                String sdf = sdfBuilder.toString();

                do
                {
                    if(line.compareTo("$$$$") == 0)
                        break;

                    if(line.compareTo(idTag) == 0)
                    {
                        line = reader.readLine();

                        id = Integer.parseInt(line.replaceAll("^" + idPrefix, ""));
                    }
                }
                while((line = reader.readLine()) != null);

                count++;
                insertStatement.setInt(1, id);
                insertStatement.setInt(2, version);
                insertStatement.setString(3, sdf);
                insertStatement.addBatch();

                if(count % batchSize == 0)
                    insertStatement.executeBatch();
            }

            if(count % batchSize != 0)
                insertStatement.executeBatch();
        }
    }


    public void loadDirectory(File directory) throws Exception
    {
        final File[] files = directory.listFiles();

        Arrays.sort(files, new Comparator<File>()
        {
            @Override
            public int compare(File arg0, File arg1)
            {
                return arg0.getName().compareTo(arg1.getName());
            }
        });;


        int version;

        try (Statement statement = connection.createStatement())
        {
            try (ResultSet rs = statement.executeQuery("select max(version)+1 from compounds"))
            {
                rs.next();
                version = rs.getInt(1);
            }
        }


        for(int i = 0; i < files.length; i++)
        {
            File file = files[i];

            if(!file.getName().endsWith(".gz"))
                continue;

            System.out.println(i + ": " + file.getName());

            InputStream fileStream = new FileInputStream(file);
            InputStream gzipStream = new GZIPInputStream(fileStream);

            parse(gzipStream, version);

            gzipStream.close();
        }

        try (PreparedStatement statement = connection.prepareStatement("delete from compounds where version < ?"))
        {
            statement.setInt(1, version);
            statement.execute();
        }
    }
}
