
#include "AtomFingerprint.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>

namespace RDKit
{
namespace AtomFingerprints
{

SparseIntVect<boost::uint32_t> *getFingerprint (const ROMol &mol,
                                                BitInfo *info)
{
	std::map<uint32_t, int> fp;
	for (auto i = mol.beginAtoms(), e = mol.endAtoms(); i != e; ++i) {
		uint32_t hsh = (*i)->getAtomicNum();
		fp[hsh] += 1;
		if (info) (*info) [hsh].insert ( (*i)->getIdx());
	}

	SparseIntVect<boost::uint32_t> *res
	    = new SparseIntVect<boost::uint32_t> (fp.empty() ?
	                                          0 : fp.rbegin()->first + 1);
	for (auto&i : fp) res->setVal (i.first, i.second);
	return res;
}

}
}
