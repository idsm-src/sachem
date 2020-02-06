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
        ArrayList<Info> sorted = new ArrayList<Info>(infos.size());

        for(SegmentCommitInfo info : infos)
            sorted.add(new Info(info, info.info.maxDoc() - context.numDeletesToMerge(info)));

        while(sorted.size() > maxSegmentCount)
        {
            sorted.sort((o1, o2) -> Long.compare(o2.docs, o1.docs));

            Info i1 = sorted.get(maxSegmentCount - 1);
            Info i2 = sorted.get(sorted.size() - 1);

            sorted.remove(sorted.size() - 1);
            i1.add(i2);
        }

        MergeSpecification specification = new MergeSpecification();

        for(Info info : sorted)
            if(info.infos.size() > 1 || info.infos.get(0).hasDeletions())
                specification.add(new OneMerge(info.infos));

        return specification;
    }


    @Override
    public MergeSpecification findMerges(MergeTrigger trigger, SegmentInfos infos, MergeContext context)
            throws IOException
    {
        return new MergeSpecification();
    }
}
