
#include "SGFingerprint.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Subgraphs/Subgraphs.h>

namespace
{

static inline void update_seed (uint32_t x, uint32_t&seed)
{
	seed ^= (uint32_t) x * 2654435761 + 2654435769
	        + (seed << 6) + (seed >> 2);
}

//unused, just for reference
#if 0
template<class I>
static inline uint32_t hashlist (const std::list<I>& l)
{
	uint32_t seed = 0;
	for (auto&&x : l)
		update_seed (x, seed);
	return seed;
}
#endif

#include <algorithm>

template<class C>
static inline uint32_t hashlist_2_sorted (uint32_t a, uint32_t b, const C& l)
{
	/*cerr << "hash: " << a << ' ' << b;
	for(auto&&i:l) cerr << ' ' << i;
	cerr << endl;*/
	std::vector<uint32_t> s (l.begin(), l.end());
	std::sort (s.begin(), s.end());
	uint32_t seed = 0;
	update_seed (a, seed);
	update_seed (b, seed);
	for (auto&i : s) update_seed (i, seed);
	return seed;
}

static inline uint32_t bond_type_hash (RDKit::Bond::BondType b)
{
	switch (b) {
	case RDKit::Bond::ZERO:
		return 0;
	case RDKit::Bond::SINGLE:
		return 1;
	case RDKit::Bond::DOUBLE:
		return 2;
	case RDKit::Bond::TRIPLE:
		return 3;
	case RDKit::Bond::QUADRUPLE:
		return 4;
	case RDKit::Bond::QUINTUPLE:
		return 5;
	case RDKit::Bond::ONEANDAHALF:
		return 6;
	case RDKit::Bond::AROMATIC:
		return 7;
	case RDKit::Bond::HYDROGEN:
		return 8;
	case RDKit::Bond::IONIC:
		return 9;
	default:
		return 10;
	}
}

static inline uint32_t atom_hash (const RDKit::Atom& a)
{
	/* TODO
	 * This might also include formal charge, but the resulting fingerprint
	 * breaks substructure matching that ignores the charge in a subset of
	 * cases.
	 */

	return a.getAtomicNum();
}

using namespace std;

using atomdesc = std::pair<uint32_t, std::list<uint32_t>>;
static bool submol_hash (const std::vector<int>&bondIds,
                         const RDKit::ROMol&mol,
                         uint32_t&result)
{
	result = 0;
	std::map<int, std::list<int>> preatoms; //atomIdx -> bondIdxs

	for (auto&bid : bondIds) {
		auto* b = mol.getBondWithIdx (bid);
		preatoms[b->getBeginAtomIdx()].push_back (bid);
		preatoms[b->getEndAtomIdx()].push_back (bid);
	}

	for (auto&i : preatoms) {
		auto*a = mol.getAtomWithIdx (i.first);
		if (a->hasProp ("dummyLabel")
		    || a->getAtomicNum() == 0
		    || a->getMass() == 0)
			return false;
	}

	//degree -> atomID -> (hash, [incoming])
	std::map<int, std::map<int, atomdesc>> atoms;

	//atomFrom -> atomTo -> hash
	std::map<int, std::map<int, uint32_t>> bonds;
	atoms[1]; //create an empty list of degree-0 atoms
	for (auto&a : preatoms) {
		auto aid = a.first;
		//pre-hash the bonds
		for (auto bid : a.second) {
			auto* b = mol.getBondWithIdx (bid);
			bonds[aid][b->getOtherAtomIdx (aid)]
			    = bond_type_hash (b->getBondType());
		}
		atoms[a.second.size()][aid]
		    = atomdesc
		      (atom_hash (*mol.getAtomWithIdx (aid)),
		       std::list<uint32_t>());
	}

	//purge the leaves until there's nothing left
	while (!atoms[1].empty()) {
		auto ats = atoms[1];
		for (auto&a : ats) {
			auto aid = a.first;
			//there's just one thing in bonds[aid]
			int addto_id = bonds[aid].begin()->first;
			uint32_t bond_hash = bonds[aid].begin()->second;
			uint32_t result_hash = hashlist_2_sorted
			                       (ats[aid].first, bond_hash,
			                        ats[aid].second);
			bonds.erase (aid);
			atoms[1].erase (aid);
			int addto_deg = bonds[addto_id].size();
			if (ats.count (addto_id)) {
				//final doublet handling!
				atomdesc other;
				other.swap (atoms[addto_deg][addto_id]);
				uint32_t other_atom_hash
				    = hashlist_2_sorted
				      (other.first, bond_hash, other.second);
				atoms[addto_deg].erase (addto_id);
				bonds[addto_id].erase (aid);
				atoms[0][addto_id].first = 666;
				atoms[0][addto_id].second.push_back
				(result_hash);
				atoms[0][addto_id].second.push_back
				(other_atom_hash);
				break;
			}
			//"normal" leaf
			atomdesc res_atom;
			res_atom.swap (atoms[addto_deg][addto_id]);
			atoms[addto_deg].erase (addto_id);
			bonds[addto_id].erase (aid);
			res_atom.second.push_back (result_hash);
			atoms[addto_deg - 1][addto_id] = res_atom;
		}
	}

	if (atoms[0].empty()) {
		//there must be a single cycle!
		for (auto&i : atoms)
			if (i.first != 2 && i.second.size())
				return false;
		//(note that the graph is connected)

		auto&ats = atoms[2];
		if (ats.empty()) return false; //this would be just weird.
		int curId = ats.begin()->first,
		    lastId = -1,
		    startId = curId;
		bool found_next = true;
		std::vector<uint32_t> cycle;
		cycle.reserve (2 * ats.size());
		atomdesc&firstAtom = ats.begin()->second;
		cycle.push_back (hashlist_2_sorted
		                 (firstAtom.first, 0, firstAtom.second));
		while (found_next) {
			found_next = false;
			for (auto&i : bonds[curId])
				if (i.first != lastId) {
					cycle.push_back (i.second);
					if (i.first == startId) break;
					atomdesc&theAtom = ats[i.first];
					cycle.push_back
					(hashlist_2_sorted
					 (theAtom.first, 0, theAtom.second));
					lastId = curId;
					curId = i.first;
					found_next = true;
					break;
				}
		}
		int minrot = 0, mindir = -1, n = cycle.size();
		for (int rot = 0; rot < n; ++rot)
			for (int dir = -1; dir <= 1; dir += 2) {
				for (int i = 0; i < n; ++i) {
					if (cycle[ (n + minrot + i * mindir) % n]
					    < cycle[ (n + rot + i * dir) % n])
						break; //must be ok
					if (cycle[ (n + minrot + i * mindir) % n]
					    == cycle[ (n + rot + i * dir) % n])
						continue; //ok
					//found better!
					minrot = rot;
					mindir = dir;
					break;
				}
			}

		uint32_t seed = 0;
		for (int i = 0; i < n; ++i)
			update_seed (cycle[ (n + minrot + i * mindir) % n], seed);

		result = seed;
		return true;
	} else {
		auto&a = *atoms[0].begin();
		result = hashlist_2_sorted
		         (a.second.first, 0, a.second.second);
		return true;
	}
}

}  // end of local anonymous namespace

namespace RDKit
{
namespace SGFingerprints
{

static void addSingleMolFP (const ROMol&mol,
                            std::map<uint32_t, int>&fp,
                            unsigned int minLen,
                            unsigned int maxLen,
                            bool forQuery,
                            BitInfo*info,
                            std::vector<int> *atomMap)
{
	/* This helps the case when multi-mol compound is disassebled with
	 * forQuery==true. Can be called on compound if forQuery is false. */

	INT_PATH_LIST_MAP allSGs; //map<int,path_list>
	if (forQuery) {
		for (int i = maxLen; i >= (int) minLen; --i) {
			allSGs = findAllSubgraphsOfLengthsMtoN
			         (mol, i, i, false);
			if (!allSGs[i].empty()) break;
		}
	} else {
		allSGs = findAllSubgraphsOfLengthsMtoN
		         (mol, minLen, maxLen, false);
	}

	PATH_LIST pl; //list<vector<int>>
	for (auto&i : allSGs) for (auto&j : i.second) {
			pl.emplace_back();
			std::swap (pl.back(), j);
		}
	allSGs.clear();

	for (auto&i : pl) {
		uint32_t hash;
		if (!submol_hash (i, mol, hash)) continue;
		fp[hash] += 1;

		if (!info) continue;
		for (auto b : i) {
			auto *bp = mol.getBondWithIdx (b);
			if (atomMap) {
				(*info) [hash].insert
				( (*atomMap) [bp->getBeginAtomIdx()]);
				(*info) [hash].insert
				( (*atomMap) [bp->getEndAtomIdx()]);
			} else {
				(*info) [hash].insert (bp->getBeginAtomIdx());
				(*info) [hash].insert (bp->getEndAtomIdx());
			}
		}
	}
}

SparseIntVect<boost::uint32_t> *getFingerprint (const ROMol &mol,
                                                unsigned int minLen,
                                                unsigned int maxLen,
                                                bool forQuery,
                                                BitInfo *info)
{
	std::map<uint32_t, int> fp;

	if (forQuery) {
		/* compound queries need to be broken to mols */
		std::vector<boost::shared_ptr<ROMol>> frags;

		std::vector<std::vector<int>> atomMaps;
		frags = RDKit::MolOps::getMolFrags
		        (mol, /*sanitize = */ false, nullptr, &atomMaps, false);

		for (size_t i = 0; i < frags.size(); ++i)
			addSingleMolFP (*frags[i], fp, minLen, maxLen,
			                forQuery, info, &atomMaps[i]);
	} else {
		addSingleMolFP (mol, fp, minLen, maxLen,
		                forQuery, info, nullptr);
	}

	SparseIntVect<boost::uint32_t> *res
	    = new SparseIntVect<boost::uint32_t> (fp.empty() ?
	                                          0 : fp.rbegin()->first + 1);
	for (auto&i : fp) res->setVal (i.first, i.second);
	return res;
}

}
}
