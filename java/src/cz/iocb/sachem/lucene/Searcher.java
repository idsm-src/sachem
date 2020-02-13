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
import java.util.concurrent.Executor;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
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
    private static ThreadFactory threadFactory = new ThreadFactory()
    {
        @Override
        public Thread newThread(Runnable r)
        {
            Thread t = Executors.defaultThreadFactory().newThread(r);
            t.setDaemon(true);
            return t;
        }
    };

    private static HashMap<String, Searcher> instances = new HashMap<String, Searcher>();

    private Path path;
    private Directory folder;
    private IndexSearcher searcher;
    private Thread watcher;
    private int threadCount;


    public static Searcher get(String name, String path, int threads) throws IOException
    {
        Searcher searcher = instances.get(name);

        if(searcher == null)
        {
            searcher = new Searcher();
            instances.put(name, searcher);
        }

        searcher.configure(path, threads);

        return searcher;
    }


    private void configure(String newPathName, int newThreadCount) throws IOException
    {
        Path newPath = Paths.get(newPathName);

        if(newThreadCount <= 1)
            newThreadCount = 0;

        if(newPath.equals(path) && threadCount == newThreadCount)
            return;

        close();


        Executor executor = newThreadCount > 1 ? Executors.newFixedThreadPool(newThreadCount, threadFactory) : null;

        threadCount = newThreadCount;
        path = newPath;
        folder = FSDirectory.open(newPath);
        searcher = new IndexSearcher(DirectoryReader.open(folder), executor);
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


    public SearchResult subsearch(byte[] molecule, int n, boolean sort, SearchMode graphMode, ChargeMode chargeMode,
            IsotopeMode isotopeMode, StereoMode stereoMode, AromaticityMode aromaticityMode, TautomerMode tautomerMode,
            int isomorphismLimit) throws IOException, CDKException, TimeoutException
    {
        Query query = new SubstructureQuery(Settings.substructureFieldName, new String(molecule), graphMode, chargeMode,
                isotopeMode, stereoMode, aromaticityMode, tautomerMode, isomorphismLimit);

        if(n == 0)
            return new SearchResult();
        else if(n >= 0)
            return new SearchResult(searcher.search(query, n));
        else if(sort)
            return new SearchResult(searcher.search(query, Integer.MAX_VALUE));
        else
            return searcher.search(query, new ResultCollectorManager());
    }


    public SearchResult simsearch(byte[] molecule, int n, boolean sort, float threshold, int depth,
            AromaticityMode aromaticityMode, TautomerMode tautomerMode)
            throws IOException, CDKException, TimeoutException
    {
        Query query = new SimilarStructureQuery(Settings.similarityFieldName, new String(molecule), threshold, depth,
                aromaticityMode, tautomerMode);


        if(n == 0)
            return new SearchResult();
        else if(n >= 0)
            return new SearchResult(searcher.search(query, n));
        else if(sort)
            return new SearchResult(searcher.search(query, Integer.MAX_VALUE));
        else
            return searcher.search(query, new ResultCollectorManager());
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
