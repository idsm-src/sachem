
#include "IOCBFingerprint.h"
#include "CRNGFingerprint.h"
#include "SGFingerprint.h"
#include "AtomFingerprint.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>

namespace
{
static inline void update_seed (const uint32_t x, uint32_t&seed)
{
	seed ^= (uint32_t) x * 2654435761 + 2654435769
	        + (seed << 6) + (seed >> 2);
}

static inline uint32_t hash3 (uint32_t a, uint32_t b, uint32_t c)
{
	uint32_t seed = 0;
	update_seed (a, seed);
	update_seed (b, seed);
	update_seed (c, seed);
	return seed;
}
}

namespace RDKit
{
namespace IOCBFingerprints
{

SparseIntVect<boost::uint32_t> *getFingerprint (const ROMol &mol,
                                                int graphSize,
                                                int circSize,
                                                int maxFeatLogCount,
                                                bool forQuery,
                                                BitInfo *info)
{
	BitInfo sgi, crngi, atomi;
	SparseIntVect<boost::uint32_t>
	*sg = SGFingerprints::getFingerprint (mol, 0, graphSize, forQuery, info ? &sgi : nullptr),
	 *crng = CRNGFingerprints::getFingerprint (mol, info ? &crngi : nullptr),
	  *atom = AtomFingerprints::getFingerprint (mol, info ? &atomi : nullptr);
	//TODO circulars

	std::map<uint32_t, int> fp;

	//TODO this is lambda.
#define process_elements(var, vari, n) \
	for (auto&i : (var)->getNonzeroElements()) { \
		uint32_t h = i.first; \
		int cnt = i.second; \
		if(forQuery) { \
			int c = 0, lc=0; \
			for (; cnt && c < maxFeatLogCount; ++c, cnt /= 2) lc=c; \
			auto hsh=hash3 ((n), lc, h); \
			fp[hsh] += 1; \
			if (info) (*info)[hsh] = vari[h]; \
		} else \
			for (int c = 0; cnt && c < maxFeatLogCount; ++c, cnt /= 2) { \
				auto hsh = hash3 ((n), c, h); \
				fp[hsh] += 1; \
				if (info) (*info)[hsh] = vari[h]; \
			} \
	}

	process_elements (sg, sgi, 1);
	process_elements (crng, crngi, 2);
	process_elements (atom, atomi, 3);

	SparseIntVect<boost::uint32_t> *res
	    = new SparseIntVect<boost::uint32_t> (fp.empty() ?
	                                          0 : fp.rbegin()->first + 1);
	for (auto&i : fp) res->setVal (i.first, i.second);

	delete sg;
	delete crng;
	delete atom;

	return res;
}

}
}
