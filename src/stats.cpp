#include <map>
#include <fstream>
#include "fingerprints/IOCBFingerprint.hpp"

extern "C"
{
#include "stats.h"
#include "fingerprints/fingerprint.h"
}


Stats *stats_create()
{
    try
    {
        return (Stats *) new std::map<uint32_t,uint32_t, std::less<uint32_t>>();
    }
    catch(...)
    {
    }

    return NULL;
}


void stats_delete(Stats *stats)
{
    try
    {
        delete stats;
    }
    catch(...)
    {
    }
}


bool stats_add(Stats *stats, const Molecule *molecule)
{
    try
    {
        std::map<uint32_t,uint32_t> &map = *((std::map<uint32_t,uint32_t> *) stats);

        std::set<uint32_t> fp = iocb_fingerprint_get(molecule, GRAPH_SIZE, RING_SIZE, MAX_FEAT_LOGCOUNT);

        for(uint32_t i : fp)
        {
            auto it = map.find(i);

            if(it == map.end())
                map[i] = 1;
            else
                map[i] = map[i] + 1;
        }

        return true;
    }
    catch(...)
    {
    }

    return false;
}


bool stats_merge(Stats *stats, StatItem *items, size_t size)
{
    try
    {
        std::map<uint32_t,uint32_t> &map = *((std::map<uint32_t,uint32_t> *) stats);

        for(size_t i = 0; i < size; i++)
        {
            auto it = map.find(items[i].fp);

            if(it == map.end())
                map[items[i].fp] = items[i].count;
            else
                map[items[i].fp] = map[items[i].fp] + items[i].count;
        }

        return true;
    }
    catch(...)
    {
    }

    return false;
}


size_t stats_get_items(Stats *stats, StatItem **items)
{
    try
    {
        std::map<uint32_t,uint32_t> &map = *((std::map<uint32_t,uint32_t> *) stats);

        *items = (StatItem *) palloc(map.size() * sizeof(StatItem));
        StatItem *item = *items;

        for(auto i : map)
        {
            item->fp = i.first;
            item->count = i.second;
            item++;
        }

        return map.size();
    }
    catch(...)
    {
    }

    *items = NULL;
    return 0;
}


bool stats_write(Stats *stats, const char *name, int limit)
{
    try
    {
        std::map<uint32_t,uint32_t> &map = *((std::map<uint32_t,uint32_t> *) stats);

        std::multimap<uint32_t,uint32_t> reverse;

        for(auto i : map)
        {
            if(i.second <= limit)
                continue;

            reverse.insert(std::make_pair(i.second, i.first));
        }

        std::ofstream stream(name, std::ios::out | std::ios::binary);

        if(stream.is_open())
        {
            for(auto i : reverse)
            {
                uint32_t value = i.second;
                stream.write((char *) &value, sizeof(uint32_t));
            }
        }

        return true;
    }
    catch(...)
    {
    }

    return false;
}
