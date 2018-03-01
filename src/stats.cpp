#include <map>
#include <fstream>
#include "fingerprints/fingerprint.h"
#include "fingerprints/IOCBFingerprint.hpp"

extern "C"
{
#include "stats.h"
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


bool stats_merge(Stats *stats, Stats *substats)
{
    try
    {
        std::map<uint32_t,uint32_t> &map = *((std::map<uint32_t,uint32_t> *) stats);
        std::map<uint32_t,uint32_t> &submap = *((std::map<uint32_t,uint32_t> *) substats);

        for(auto i : submap)
        {
            auto it = map.find(i.first);

            if(it == map.end())
                map[i.first] = i.second;
            else
                map[i.first] = map[i.first] + i.second;
        }

        return true;
    }
    catch(...)
    {
    }

    return false;
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
