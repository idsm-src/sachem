package cz.iocb.sachem.lucene;

import static java.nio.file.StandardWatchEventKinds.ENTRY_DELETE;
import static java.nio.file.StandardWatchEventKinds.OVERFLOW;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.util.HashMap;
import java.util.concurrent.TimeoutException;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.similarities.BooleanSimilarity;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;
import org.openscience.cdk.exception.CDKException;
import cz.iocb.sachem.molecule.AromaticityMode;
import cz.iocb.sachem.molecule.ChargeMode;
import cz.iocb.sachem.molecule.IsotopeMode;
import cz.iocb.sachem.molecule.SearchMode;
import cz.iocb.sachem.molecule.StereoMode;
import cz.iocb.sachem.molecule.TautomerMode;



public class Searcher
{
    private static HashMap<String, Searcher> instances = new HashMap<String, Searcher>();

    private Path path;
    private Directory folder;
    private IndexSearcher searcher;
    private Thread watcher;


    public static Searcher get(String name)
    {
        Searcher searcher = instances.get(name);

        if(searcher == null)
        {
            searcher = new Searcher();
            instances.put(name, searcher);
        }

        return searcher;
    }


    public void setFolder(String newPathName) throws IOException
    {
        Path newPath = Paths.get(newPathName);

        if(newPath.equals(path))
            return;

        close();

        path = newPath;
        folder = FSDirectory.open(newPath);
        searcher = new IndexSearcher(DirectoryReader.open(folder));
        searcher.setSimilarity(new BooleanSimilarity());


        watcher = new Thread()
        {
            @Override
            public void run()
            {
                try
                {
                    IndexSearcher searcher = Searcher.this.searcher;

                    WatchService service = FileSystems.getDefault().newWatchService();
                    path.getParent().register(service, ENTRY_DELETE);

                    while(true)
                    {
                        WatchKey key = service.take();

                        for(WatchEvent<?> event : key.pollEvents())
                        {
                            if(event.kind() != OVERFLOW && event.context().equals(path.getFileName()))
                            {
                                //TODO: add check whether it is used

                                if(searcher == Searcher.this.searcher)
                                    close();

                                return;
                            }
                        }

                        if(!key.reset())
                            break;
                    }
                }
                catch(IOException e)
                {
                    e.printStackTrace();
                }
                catch(InterruptedException e)
                {
                    return;
                }
            }
        };

        watcher.setDaemon(true);
        watcher.start();
    }


    public int indexSize()
    {
        return searcher.getIndexReader().numDocs();
    }


    public SearcherHandler subsearch(byte[] molecule, SearchMode graphMode, ChargeMode chargeMode,
            IsotopeMode isotopeMode, StereoMode stereoMode, AromaticityMode aromaticityMode, TautomerMode tautomerMode,
            int bufferSize) throws IOException, CDKException, TimeoutException
    {
        Query query = new SubstructureQuery(Settings.substructureFieldName, new String(molecule), graphMode, chargeMode,
                isotopeMode, stereoMode, aromaticityMode, tautomerMode);

        return new SearcherHandler(searcher, query, bufferSize);
    }


    public SearcherHandler simsearch(byte[] molecule, float threshold, int depth, AromaticityMode aromaticityMode,
            TautomerMode tautomerMode, int bufferSize) throws IOException, CDKException, TimeoutException
    {
        Query query = new SimilarStructureQuery(Settings.similarityFieldName, new String(molecule), threshold, depth,
                aromaticityMode, tautomerMode);

        return new SearcherHandler(searcher, query, bufferSize);
    }


    protected void close() throws IOException
    {
        path = null;

        try
        {
            if(searcher != null)
                searcher.getIndexReader().close();
        }
        finally
        {
            searcher = null;
        }


        try
        {
            if(watcher != null)
                watcher.interrupt();
        }
        finally
        {
            watcher = null;
        }


        try
        {
            if(folder != null)
                folder.close();
        }
        finally
        {
            folder = null;
        }
    }
}
