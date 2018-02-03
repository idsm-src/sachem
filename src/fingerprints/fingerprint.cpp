#include <map>
#include <vector>
#include <fstream>
#include "IOCBFingerprint.hpp"

extern "C"
{
#include "molecule.h"
#include "fingerprint.h"

PG_MODULE_MAGIC;
}


#define GRAPH_SIZE              7
#define RING_SIZE               3
#define MAX_FEAT_LOGCOUNT       5
#define QUERY_ATOM_COVERAGE     2
#define QUERY_MAX_FPS          32


static const unsigned char b64str[65] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
static std::map<uint32_t, int> fporder __attribute__ ((init_priority(200)));


void __attribute__ ((constructor(300))) fingerprint_init(void)
{
    std::ifstream stream(DATADIR "/fporder.bin", std::ios::in | std::ios::binary);

    if(stream.is_open())
    {
        int fpn = 0;

        while(!stream.eof())
        {
            uint32_t fp;
            stream.read((char *) &fp, sizeof(uint32_t));
            fp = ntohl(fp);

            fporder[fp] = fpn++;
        }
    }
}


static inline void write_bitword(char *buffer, uint32_t fp)
{
    for(int i = 0; i < 6; i++)
        buffer[i] = b64str[(fp >> (6 * i)) & 0x3f];
}


Fingerprint fingerprint_get(const Molecule *molecule, void *(*alloc)(size_t))
{
    try
    {
        std::set<uint32_t> res = iocb_fingerprint_get(molecule, GRAPH_SIZE, RING_SIZE, MAX_FEAT_LOGCOUNT);

        if(res.size() > 0)
        {
            size_t size = res.size() * 7;
            char *data = (char *) (*alloc)(size);

            if(data == NULL)
                return {.size = (size_t) -1, .data = NULL};


            Fingerprint fp = {size : size - 1, data: data};

            for(uint32_t i : res)
            {
                write_bitword(data, i);
                data[6] = ' ';
                data += 7;
            }

            data[-1] = '\0';
            return fp;
        }
        else
        {
            return {.size = 0, .data = NULL};
        }
    }
    catch(...)
    {
    }

    return {.size = (size_t) -1, .data = NULL};
}


Fingerprint fingerprint_get_query(const Molecule *molecule, void *(*alloc)(size_t))
{
    try
    {
        BitInfo info;
        std::set<uint32_t> res = iocb_fingerprint_get(molecule, GRAPH_SIZE, RING_SIZE, MAX_FEAT_LOGCOUNT, true, &info);


        // convert and pre-sort the fingerprints
        std::map<int, uint32_t> fpi;
        int unknownId = -1;

        for(uint32_t i : res)
        {
            auto o = fporder.find(i);

            if(o == fporder.end())
                // if the fingerprint is not known to fporder (which it should be but keeping that database in shape
                // is not very easy), let's assume it's very good (and put it on the beginning of the queue...
                fpi[unknownId--] = i;
            else
                fpi[o->second] = i;
        }


        // find a decent coverage
        std::set<uint32_t> fps;

        std::vector<int> coverage;
        int uncovered = molecule->atomCount;
        int nfps = 0;
        coverage.resize(uncovered, 0);

        for(auto &i : fpi)
        {
            if(uncovered <= 0)
                break;

            if(nfps >= QUERY_MAX_FPS)
                break;

            bool found = false;

            for(auto a : info[i.second])
            {
                if(coverage[a] < QUERY_ATOM_COVERAGE)
                {
                    found = true;
                    coverage[a]++;

                    if(coverage[a] == QUERY_ATOM_COVERAGE)
                        uncovered--;
                }
            }

            if(found)
            {
                fps.insert(i.second);
                nfps++;
            }
        }


        if(fps.size() > 0)
        {
            size_t size = fps.size() * 9;
            char *data = (char *) (*alloc)(size);

            if(data == NULL)
                return {size : (size_t) -1, data : NULL};

            Fingerprint fp = {size : size - 1, data: data};

            for(uint32_t i : fps)
            {
                data[0] = '"';
                write_bitword(data + 1, i);
                data[7] = '"';
                data[8] = ' ';
                data += 9;
            }

            data[-1] = '\0';
            return fp;
        }
        else
        {
            return {.size = 0, .data = NULL};
        }
    }
    catch(...)
    {
    }

    return {.size = (size_t) -1, .data = NULL};
}
