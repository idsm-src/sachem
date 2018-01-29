#include "IOCBFingerprint.hpp"
#include "AtomFingerprint.hpp"
#include "CRNGFingerprint.hpp"
#include "SGFingerprint.hpp"


static inline void update_seed (const uint32_t x, uint32_t&seed)
{
	seed ^= (uint32_t) x * 2654435761 + 2654435769 + (seed << 6) + (seed >> 2);
}


static inline uint32_t hash3 (uint32_t a, uint32_t b, uint32_t c)
{
	uint32_t seed = 0;
	update_seed (a, seed);
	update_seed (b, seed);
	update_seed (c, seed);
	return seed;
}


static inline void process_elements(std::map<uint32_t, int> &var, BitInfo &vari, int n, std::set<uint32_t> &fp,
        int maxFeatLogCount, bool forQuery, BitInfo *info)
{
    for(auto&i : (var))
    {
        uint32_t h = i.first;
        int cnt = i.second;

        if(forQuery)
        {
            int lc=0;

            for(int c = 0; cnt && c < maxFeatLogCount; ++c, cnt /= 2)
                lc=c;

            auto hsh = hash3(n, lc, h);
            fp.insert(hsh);

            if(info)
                (*info)[hsh] = vari[h];
        }
        else
        {
            for(int c = 0; cnt && c < maxFeatLogCount; ++c, cnt /= 2)
            {
                auto hsh = hash3(n, c, h);
                fp.insert(hsh);

                if(info)
                    (*info)[hsh] = vari[h];
            }
        }
    }
}


std::set<uint32_t> iocb_fingerprint_get(const Molecule *molecule, int graphSize, int circSize, int maxFeatLogCount,
        bool forQuery, BitInfo *info)
{
    std::set<uint32_t> fp;

    BitInfo sgi;
    std::map<uint32_t, int> sg = sg_fingerprint_get(molecule, 0, graphSize, forQuery, info ? &sgi : nullptr);
    process_elements (sg, sgi, 1, fp, maxFeatLogCount, forQuery, info);

	BitInfo crngi;
	std::map<uint32_t, int> crng = crng_fingerprint_get(molecule, info ? &crngi : nullptr);
    process_elements (crng, crngi, 2, fp, maxFeatLogCount, forQuery, info);

	BitInfo atomi;
	std::map<uint32_t, int> atom = atom_fingerprint_get(molecule, info ? &atomi : nullptr);
	process_elements (atom, atomi, 3, fp, maxFeatLogCount, forQuery, info);

	return fp;
}
