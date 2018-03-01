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
            fporder[fp] = fpn++;
        }
    }
}


static inline void write_bitword(char *buffer, uint32_t fp)
{
    for(int i = 0; i < 6; i++)
        buffer[i] = b64str[(fp >> (6 * i)) & 0x3f];
}


static inline IntegerFingerprint integer_fingerprint_create(std::set<uint32_t> &res, void *(*alloc)(size_t))
{
    if(res.size() > 0)
    {
        size_t size = res.size();
        int32_t *data = (int *) (*alloc)(size * sizeof(int32_t));

        if(data == NULL)
            return {.size = (size_t) -1, .data = NULL};


        IntegerFingerprint fp = {size : size, data: data};

        for(uint32_t i : res)
            *(data++) = i;

        return fp;
    }
    else
    {
        return {.size = 0, .data = NULL};
    }
}


static inline std::set<uint32_t> fingerprint_get_native(const Molecule *molecule)
{
    return iocb_fingerprint_get(molecule, GRAPH_SIZE, RING_SIZE, MAX_FEAT_LOGCOUNT);
}


static inline std::set<uint32_t> fingerprint_get_query_native(const Molecule *molecule)
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

    return fps;
}


StringFingerprint string_fingerprint_get(const Molecule *molecule, void *(*alloc)(size_t))
{
    try
    {
        std::set<uint32_t> res = fingerprint_get_native(molecule);

        if(res.size() > 0)
        {
            size_t size = res.size() * 7;
            char *data = (char *) (*alloc)(size);

            if(data == NULL)
                return {.size = (size_t) -1, .data = NULL};


            StringFingerprint fp = {size : size - 1, data: data};

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


StringFingerprint string_fingerprint_get_query(const Molecule *molecule, void *(*alloc)(size_t))
{
    try
    {
        std::set<uint32_t> fps = fingerprint_get_query_native(molecule);

        if(fps.size() > 0)
        {
            size_t size = fps.size() * 9;
            char *data = (char *) (*alloc)(size);

            if(data == NULL)
                return {size : (size_t) -1, data : NULL};

            StringFingerprint fp = {size : size - 1, data: data};

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


IntegerFingerprint integer_fingerprint_get(const Molecule *molecule, void *(*alloc)(size_t))
{
    try
    {
        std::set<uint32_t> res = fingerprint_get_native(molecule);
        return integer_fingerprint_create(res, alloc);
    }
    catch(...)
    {
    }

    return {.size = (size_t) -1, .data = NULL};
}


IntegerFingerprint integer_fingerprint_get_query(const Molecule *molecule, void *(*alloc)(size_t))
{
    try
    {
        std::set<uint32_t> res = fingerprint_get_query_native(molecule);
        return integer_fingerprint_create(res, alloc);
    }
    catch(...)
    {
    }

    return {.size = (size_t) -1, .data = NULL};
}
