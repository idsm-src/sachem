package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.nio.file.Paths;
import org.apache.lucene.index.SegmentCommitInfo;
import org.apache.lucene.index.SegmentInfos;
import org.apache.lucene.store.FSDirectory;



public class IndexInfo
{
    static class SegmentInfo
    {
        public String name;
        public int molecules;
        public int deletes;
        public long size;

        public SegmentInfo(SegmentCommitInfo segment) throws IOException
        {
            this.name = segment.info.name;
            this.molecules = segment.info.maxDoc() - segment.getDelCount();
            this.deletes = segment.getDelCount();
            this.size = segment.sizeInBytes();
        }
    }


    public static SegmentInfo[] getSegmentInfos(String path) throws IOException
    {
        FSDirectory folder = FSDirectory.open(Paths.get(path));

        SegmentInfos infos = SegmentInfos.readLatestCommit(folder);
        SegmentInfo[] result = new SegmentInfo[infos.size()];

        int i = 0;

        for(SegmentCommitInfo info : infos)
            result[i++] = new SegmentInfo(info);

        return result;
    }
}
