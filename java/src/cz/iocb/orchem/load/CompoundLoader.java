package cz.iocb.orchem.load;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.util.Arrays;
import java.util.Comparator;
import java.util.zip.GZIPInputStream;



public class CompoundLoader
{
    private static final int batchSize = 1000;


    private static void parse(InputStream inputStream, String idTag, String idPrefix) throws Exception
    {
        Reader decoder = new InputStreamReader(inputStream, "US-ASCII");
        BufferedReader reader = new BufferedReader(decoder);

        String line;

        try (Connection connection = ConnectionPool.getConnection())
        {
            try (PreparedStatement insertStatement = connection
                    .prepareStatement("insert into compounds (id, molfile) values (?,?)"))
            {
                int count = 0;

                while((line = reader.readLine()) != null)
                {
                    int id = -1;
                    StringBuilder sdfBuilder = new StringBuilder();


                    sdfBuilder.append('\n');
                    sdfBuilder.append(line);

                    line = reader.readLine();
                    sdfBuilder.append(line);
                    sdfBuilder.append('\n');



                    for(int i = 0; i < 2; i++)
                    {
                        line = reader.readLine();
                        sdfBuilder.append(line);
                        sdfBuilder.append('\n');
                    }

                    int atoms = Integer.parseInt(line.substring(0, 3).trim());
                    int bonds = Integer.parseInt(line.substring(3, 6).trim());


                    for(int i = 0; i < atoms; i++)
                    {
                        line = reader.readLine();
                        sdfBuilder.append(line);
                        sdfBuilder.append('\n');
                    }


                    for(int i = 0; i < bonds; i++)
                    {
                        line = reader.readLine();
                        sdfBuilder.append(line);
                        sdfBuilder.append('\n');
                    }

                    while((line = reader.readLine()) != null)
                    {
                        sdfBuilder.append(line);
                        sdfBuilder.append('\n');

                        if(line.compareTo("M  END") == 0)
                            break;
                    }

                    String sdf = sdfBuilder.toString();

                    while((line = reader.readLine()) != null)
                    {
                        if(line.compareTo("$$$$") == 0)
                            break;

                        if(line.compareTo(idTag) == 0)
                        {
                            line = reader.readLine();

                            id = Integer.parseInt(line.replaceAll("^" + idPrefix, ""));
                        }
                    }

                    count++;
                    insertStatement.setInt(1, id);
                    insertStatement.setString(2, sdf);
                    insertStatement.addBatch();

                    if(count % batchSize == 0)
                        insertStatement.executeBatch();
                }

                if(count % batchSize != 0)
                    insertStatement.executeBatch();
            }
        }
    }


    public static void loadDirectory(File directory, String idTag, String idPrefix) throws Exception
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


        for(int i = 0; i < files.length; i++)
        {
            File file = files[i];

            if(!file.getName().endsWith(".gz"))
                continue;

            System.out.println(i + ": " + file.getName());

            InputStream fileStream = new FileInputStream(file);
            InputStream gzipStream = new GZIPInputStream(fileStream);

            parse(gzipStream, idTag, idPrefix);

            gzipStream.close();
        }
    }
}
