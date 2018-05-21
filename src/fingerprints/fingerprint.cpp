#include <map>
#include <vector>
#include <fstream>
#include "IOCBFingerprint.hpp"

extern "C"
{
#include "fporder.h"
#include "molecule.h"
#include "fingerprint.h"
#include "sachem.h"
}


#define QUERY_ATOM_COVERAGE     2
#define QUERY_MAX_FPS          32


static bool initialized = false;
static const unsigned char b64str[65] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
static std::map<uint32_t, int> fporder;


void fingerprint_init(void)
{
    char *fporderPath = get_file_path(ORDER_FILE);
    std::ifstream stream(fporderPath, std::ios::in | std::ios::binary);

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
    else
    {
        elog(WARNING, "cannot load the fporder file: %s", fporderPath);
    }
}


static inline void write_bitword(char *buffer, uint32_t fp)
{
    for(int i = 0; i < 6; i++)
        buffer[i] = b64str[(fp >> (6 * i)) & 0x3f];
}


static inline IntegerFingerprint integer_fingerprint_create(std::set<uint32_t> &res)
{
    if(res.size() > 0)
    {
        size_t size = res.size();
        int32_t *data = (int *) palloc_extended(size * sizeof(int32_t), MCXT_ALLOC_NO_OOM);

        if(data == NULL)
            throw std::bad_alloc();


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


static inline std::set<uint32_t> substructure_fingerprint_get_native(const Molecule *molecule)
{
    return iocb_substructure_fingerprint_get(molecule, GRAPH_SIZE, MAX_FEAT_LOGCOUNT);
}


static inline std::set<uint32_t> substructure_fingerprint_get_query_native(const Molecule *molecule)
{
    if(unlikely(initialized == false))
    {
        fingerprint_init();
        initialized = true;
    }


    BitInfo info;
    std::set<uint32_t> res = iocb_substructure_fingerprint_get(molecule, GRAPH_SIZE, MAX_FEAT_LOGCOUNT, true, &info);


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


static inline std::set<uint32_t> similarity_fingerprint_get_native(const Molecule *molecule)
{
    return iocb_similarity_fingerprint_get(molecule, CIRC_SIZE, MAX_FEAT_LOGCOUNT);
}


StringFingerprint string_substructure_fingerprint_get(const Molecule *molecule)
{
    SAFE_CPP_BEGIN;

    std::set<uint32_t> res = substructure_fingerprint_get_native(molecule);

    if(res.size() > 0)
    {
        size_t size = res.size() * 7;
        char *data = (char *) palloc_extended(size, MCXT_ALLOC_NO_OOM);

        if(data == NULL)
            throw std::bad_alloc();


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

    SAFE_CPP_END;
}


StringFingerprint string_substructure_fingerprint_get_query(const Molecule *molecule)
{
    SAFE_CPP_BEGIN;

    std::set<uint32_t> fps = substructure_fingerprint_get_query_native(molecule);

    if(fps.size() > 0)
    {
        size_t size = fps.size() * 9;
        char *data = (char *) palloc_extended(size, MCXT_ALLOC_NO_OOM);

        if(data == NULL)
            throw std::bad_alloc();

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

    SAFE_CPP_END;
}


IntegerFingerprint integer_substructure_fingerprint_get(const Molecule *molecule)
{
    SAFE_CPP_BEGIN;

    std::set<uint32_t> res = substructure_fingerprint_get_native(molecule);
    return integer_fingerprint_create(res);

    SAFE_CPP_END;
}


IntegerFingerprint integer_substructure_fingerprint_get_query(const Molecule *molecule)
{
    SAFE_CPP_BEGIN;

    std::set<uint32_t> res = substructure_fingerprint_get_query_native(molecule);
    return integer_fingerprint_create(res);

    SAFE_CPP_END;
}


IntegerFingerprint integer_similarity_fingerprint_get(const Molecule *molecule)
{
    SAFE_CPP_BEGIN;

    std::set<uint32_t> res = similarity_fingerprint_get_native(molecule);
    return integer_fingerprint_create(res);

    SAFE_CPP_END;
}


IntegerFingerprint integer_similarity_fingerprint_get_query(const Molecule *molecule)
{
    SAFE_CPP_BEGIN;

    std::set<uint32_t> res = similarity_fingerprint_get_native(molecule);
    return integer_fingerprint_create(res);

    SAFE_CPP_END;
}
