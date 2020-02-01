package cz.iocb.sachem.lucene;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import org.apache.lucene.index.MergePolicy;
import org.apache.lucene.index.MergeTrigger;
import org.apache.lucene.index.SegmentCommitInfo;
import org.apache.lucene.index.SegmentInfos;



public class BalancedMergePolicy extends MergePolicy
{
    private static class Info
    {
        List<SegmentCommitInfo> infos = new LinkedList<SegmentCommitInfo>();
        long docs = 0;

        public Info(SegmentCommitInfo info, int docs)
        {
            this.infos.add(info);
            this.docs = docs;
        }

        void add(Info o)
        {
            infos.addAll(o.infos);
            docs += o.docs;
        }
    }


    private final int segmentsCount;


    public BalancedMergePolicy(int segmentsCount)
    {
        this.segmentsCount = segmentsCount;
    }


    @Override
    public MergeSpecification findForcedDeletesMerges(SegmentInfos infos, MergeContext context) throws IOException
    {
        MergeSpecification specification = new MergeSpecification();

        for(SegmentCommitInfo info : infos)
            if(info.hasDeletions())
                specification.add(new OneMerge(Collections.singletonList(info)));

        return specification;
    }


    @Override
    public MergeSpecification findForcedMerges(SegmentInfos infos, int maxSegmentCount,
            Map<SegmentCommitInfo, Boolean> segmentsToMerge, MergeContext context) throws IOException
    {
        throw new UnsupportedOperationException();
    }


    @Override
    public MergeSpecification findMerges(MergeTrigger trigger, SegmentInfos infos, MergeContext context)
            throws IOException
    {
        ArrayList<Info> sorted = new ArrayList<Info>(infos.size());

        for(SegmentCommitInfo info : infos)
            sorted.add(new Info(info, info.info.maxDoc() - context.numDeletesToMerge(info)));

        while(sorted.size() > segmentsCount)
        {
            sorted.sort((o1, o2) -> Long.compare(o2.docs, o1.docs));

            Info i1 = sorted.get(segmentsCount - 1);
            Info i2 = sorted.get(sorted.size() - 1);

            sorted.remove(sorted.size() - 1);
            i1.add(i2);
        }

        MergeSpecification specification = new MergeSpecification();

        for(Info info : sorted)
            if(info.infos.size() > 1)
                specification.add(new OneMerge(info.infos));

        return specification;
    }
}
