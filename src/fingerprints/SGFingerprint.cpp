/*
 * Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
 * Copyright (C) 2016-2017 Miroslav Kratochvil
 * Copyright (C) 2016-2017 Jakub Galgonek
 *
 * Subgraph functions of this file are part of the RDKit.
 *
 * The contents are covered by the terms of the BSD license which is included in the file license.txt, found at the root
 * of the RDKit source tree.
 */
#include <list>
#include <vector>
#include <algorithm>
#include "SGFingerprint.hpp"

extern "C"
{
#include "molecule.h"
}


typedef std::vector<int> PathType;
typedef std::list<PathType> PathList;
typedef std::map<uint32_t, std::set<uint32_t>> BitInfo;
typedef std::pair<uint32_t, std::list<uint32_t>> AtomDesc;


static inline void update_seed(uint32_t x, uint32_t &seed)
{
	seed ^= (uint32_t) x * 2654435761 + 2654435769 + (seed << 6) + (seed >> 2);
}


template<class C>
static inline uint32_t list_hash(uint32_t a, uint32_t b, const C &l)
{
	std::vector<uint32_t> s(l.begin(), l.end());
	std::sort(s.begin(), s.end());

	uint32_t seed = 0;
	update_seed(a, seed);
	update_seed(b, seed);

	for(uint32_t &i : s)
	    update_seed(i, seed);

	return seed;
}


static inline uint32_t bond_hash(const Molecule *molecule, int b)
{
    return molecule_get_bond_type(molecule, b);
}


static inline uint32_t atom_hash(const Molecule *molecule, int a)
{
	return molecule_get_atom_number(molecule, a);
}


static bool submol_hash(const std::vector<int> &bondIds, const Molecule *molecule, uint32_t &result)
{
    result = 0;

	std::map<int, std::list<int>> preatoms;         // atomIdx -> bondIdxs

	for(const int bid : bondIds)
	{
        preatoms[molecule_bond_atoms(molecule, bid)[0]].push_back(bid);
        preatoms[molecule_bond_atoms(molecule, bid)[1]].push_back(bid);
    }

	std::map<int, std::map<int, AtomDesc>> atoms;   // degree -> atomID -> (hash, [incoming])
	std::map<int, std::map<int, uint32_t>> bonds;   // atomFrom -> atomTo -> hash

	atoms[1]; // create an empty list of degree-0 atoms

	for(auto &a : preatoms)
	{
		int aid = a.first;

		// pre-hash the bonds
		for(int bid : a.second)
		    bonds[aid][molecule_get_bond_connected_atom(molecule, bid, aid)] = bond_hash(molecule, bid);

        atoms[a.second.size()][aid] = AtomDesc(atom_hash(molecule, aid), std::list<uint32_t>());
	}

	// purge the leaves until there is nothing left
	while(!atoms[1].empty())
	{
	    std::map<int, AtomDesc> ats = atoms[1];

		for(auto &a : ats)
		{
			int aid = a.first;

			// there is just one thing in bonds[aid]
			int addto_id = bonds[aid].begin()->first;
			uint32_t bond_hash = bonds[aid].begin()->second;

			uint32_t result_hash = list_hash(ats[aid].first, bond_hash, ats[aid].second);

			bonds.erase(aid);
			atoms[1].erase(aid);
			int addto_deg = bonds[addto_id].size();

			if(ats.count(addto_id))
			{
				// final doublet handling!
				AtomDesc other;
				other.swap(atoms[addto_deg][addto_id]);

				uint32_t other_atom_hash = list_hash(other.first, bond_hash, other.second);
				atoms[addto_deg].erase(addto_id);
				bonds[addto_id].erase(aid);
				atoms[0][addto_id].first = 666; //FIXME: use different value
				atoms[0][addto_id].second.push_back(result_hash);
				atoms[0][addto_id].second.push_back(other_atom_hash);
				break;
			}

			// "normal" leaf
			AtomDesc res_atom;
			res_atom.swap(atoms[addto_deg][addto_id]);
			atoms[addto_deg].erase(addto_id);
			bonds[addto_id].erase(aid);
			res_atom.second.push_back(result_hash);
			atoms[addto_deg - 1][addto_id] = res_atom;
		}
	}

	if(atoms[0].empty())
	{
		// there must be a single cycle!
		for(auto &i : atoms)
			if(i.first != 2 && i.second.size())
				return false;
		// (note that the graph is connected)

		auto &ats = atoms[2];

		if(ats.empty())
		    return false; // this would be just weird.

		int curId = ats.begin()->first;
		int lastId = -1;
		int startId = curId;
		bool found_next = true;
		std::vector<uint32_t> cycle;

		cycle.reserve(2 * ats.size());
		AtomDesc &firstAtom = ats.begin()->second;
		cycle.push_back(list_hash(firstAtom.first, 0, firstAtom.second));

		while(found_next)
		{
			found_next = false;

			for(auto &i : bonds[curId])
			{
				if(i.first != lastId)
				{
					cycle.push_back(i.second);

					if(i.first == startId)
					    break;

					AtomDesc &theAtom = ats[i.first];
					cycle.push_back(list_hash(theAtom.first, 0, theAtom.second));
					lastId = curId;
					curId = i.first;
					found_next = true;
					break;
				}
			}
		}

		int minrot = 0;
		int mindir = -1;
		int n = cycle.size();

		for(int rot = 0; rot < n; rot++)
		{
			for(int dir = -1; dir <= 1; dir += 2)
			{
				for(int i = 0; i < n; i++)
				{
					if(cycle[(n + minrot + i * mindir) % n] < cycle[(n + rot + i * dir) % n])
						break; // must be ok

					if(cycle[(n + minrot + i * mindir) % n] == cycle[(n + rot + i * dir) % n])
						continue; // ok

					// found better!
					minrot = rot;
					mindir = dir;
					break;
				}
			}
		}

		uint32_t seed = 0;

		for(int i = 0; i < n; i++)
			update_seed(cycle[(n + minrot + i * mindir) % n], seed);

		result = seed;
		return true;
	}
	else
	{
		auto &a = *atoms[0].begin();
		result = list_hash(a.second.first, 0, a.second.second);
		return true;
	}
}


static void get_neighbor_list(const Molecule *molecule, std::map<int, std::vector<int>> &nbrs)
{
    nbrs.clear();

    int nAtoms = molecule->atomCount;

    // create a list of neighbors for each bond
    for(int i = 0; i < nAtoms; i++)
    {
        if(molecule_get_atom_number(molecule, i) <= H_ATOM_NUMBER)
            continue;

        int size = molecule_get_bond_list_size(molecule, i);
        int *bondedAtoms = molecule_get_bond_list(molecule, i);

        for(int l = 0; l < size; l++)
        {
            if(molecule_get_atom_number(molecule, bondedAtoms[l]) <= H_ATOM_NUMBER)
                continue;

            int bid1 = molecule_get_bond(molecule, i, bondedAtoms[l]);

            if(nbrs.find(bid1) == nbrs.end())
            {
                std::vector<int> nlst;
                nbrs[bid1] = nlst;
            }

            for(int k = 0; k < size; k++)
            {
                if(molecule_get_atom_number(molecule, bondedAtoms[k]) <= H_ATOM_NUMBER)
                    continue;

                int bid2 = molecule_get_bond(molecule, i, bondedAtoms[k]);

                if(bid1 != bid2)
                    nbrs[bid1].push_back(bid2);  //FIXME: pathListType should probably be container of pointers?
            }
        }
    }
}


/**
 * @param nbrs      neighbors for each bond
 * @param spath     the current path to be build upon
 * @param cands     neighbors of current path
 * @param lowerLen  lower limit of the subgraph lengths we are interested in
 * @param upperLen  the maximum subgraph len we are interested in
 * @param forbidden bonds that have been covered already we do not want reference passing for forbidden
 *                  it gets altered through the processand we want fresh start everytime we buble back
 *                  up to "FindAllSubGraphs"
 * @param res       the final list of subgraphs
 */
static void recurse_walk_range(std::map<int, std::vector<int>> &nbrs, PathType &spath, std::vector<int> &cands,
        uint lowerLen, uint upperLen, std::vector<uint8_t> forbidden, std::map<int, PathList> &res)
{
    uint nsize = spath.size();

    if(nsize >= lowerLen && nsize <= upperLen)
        res[nsize].push_back(spath);

    // end case for recursion
    if(nsize == upperLen)
        return;


    // if the path is already bigger than desired size
    if(nsize > upperLen)
        return;


    // we  have the candidates that can be used to add to the existing path try extending the subgraphs
    while(cands.size() != 0)
    {
        int next = cands.back();  // start with the last one in the candidate list
        cands.pop_back();

        if(!forbidden[next])
        {
            // this bond should not appear in the later subgraphs
            forbidden[next] = 1;

            // update a local stack before the next recursive call
            std::vector<int> tstack = cands;

            for(std::vector<int>::iterator bid = nbrs[next].begin(); bid != nbrs[next].end(); bid++)
                if(!forbidden[*bid])
                    tstack.push_back(*bid);

            PathType tpath = spath;
            tpath.push_back(next);

            recurse_walk_range(nbrs, tpath, tstack, lowerLen, upperLen, forbidden, res);
        }
    }
}


static std::map<int, PathList> find_subgraphs(const Molecule *molecule, uint lowerLen, uint upperLen, int rootedAtAtom)
{
    std::vector<uint8_t> forbidden(molecule->bondCount);

    std::map<int, std::vector<int>> nbrs;
    get_neighbor_list(molecule, nbrs);

    // start path at each bond
    std::map<int, PathList> res;

    for(uint idx = lowerLen; idx <= upperLen; idx++)
    {
        PathList ordern;
        res[idx] = ordern;
    }

    // start paths at each bond:
    for(std::map<int, std::vector<int>>::iterator nbi = nbrs.begin(); nbi != nbrs.end(); nbi++)
    {
        int i = (*nbi).first;

        // if we are only returning paths rooted at a particular atom, check now that this bond involves that atom:
        if(rootedAtAtom >= 0 && molecule_bond_atoms(molecule, i)[0] != rootedAtAtom &&
                molecule_bond_atoms(molecule, i)[1] != rootedAtAtom)
            continue;

        // do not come back to this bond in the later subgraphs
        if(forbidden[i])
            continue;

        forbidden[i] = 1;

        // start the recursive path building with the current bond
        PathType spath;
        spath.clear();
        spath.push_back(i);

        // neighbors of this bond are the next candidates
        std::vector<int> cands = nbrs[i];

        // now call the recursive function little bit different from the python version
        // the result list of paths is passed as a reference, instead of on the fly appending
        recurse_walk_range(nbrs, spath, cands, lowerLen, upperLen, forbidden, res);
    }

    nbrs.clear();

    return res;  //FIXME: need some verbose testing code here
}


static void add_molecule_fp(const Molecule *molecule, std::map<uint32_t, int> &fp, uint minLen, uint maxLen, BitInfo*info)
{
	std::map<int, PathList> allSGs = find_subgraphs(molecule, minLen, maxLen, -1);


    for(auto &i : allSGs)
    {
        for(auto &j : i.second)
        {
            uint32_t hash;

            if(!submol_hash(j, molecule, hash))
                continue;

            fp[hash] += 1;


            if(!info)
                continue;

            for(auto b : j)
            {
                (*info)[hash].insert(molecule_bond_atoms(molecule, b)[0]);
                (*info)[hash].insert(molecule_bond_atoms(molecule, b)[1]);
            }
        }
	}
}


static inline void fragment_walk(const Molecule *molecule, int atom, std::vector<int> &visitedAtoms,
        std::vector<int> &visitedBonds, int &visited)
{
    visitedAtoms[atom] = 1;

    int size = molecule_get_bond_list_size(molecule, atom);
    int *neighbors = molecule_get_bond_list(molecule, atom);

    for(int i = 0; i < size; i++)
    {
        int other = neighbors[i];

        int bond = molecule_get_bond(molecule, atom, other);

        if(!visitedBonds[bond])
        {
            visitedBonds[bond] = 1;
            visited++;
        }

        if(!visitedAtoms[other])
            fragment_walk(molecule, other, visitedAtoms, visitedBonds, visited);
    }
}


std::map<uint32_t, int> sg_fingerprint_get(const Molecule *molecule, uint minLen, uint maxLen, bool forQuery, BitInfo *info)
{
	if(forQuery)
	{
	    int minQueryLen = maxLen;

	    std::vector<int> visitedAtoms(molecule->atomCount);
	    std::vector<int> visitedBonds(molecule->bondCount);

	    for(int a = 0; a < molecule->atomCount; a++)
	    {
	        if(!visitedAtoms[a])
	        {
	            int visited = 0;
	            fragment_walk(molecule, a, visitedAtoms, visitedBonds, visited);

	            if(visited > 0 && visited < minQueryLen && visited >= minLen)
	                minQueryLen = visited;
	        }
	    }

	    minLen = minQueryLen;
	}


    std::map<uint32_t, int> fp;
    add_molecule_fp(molecule, fp, minLen, maxLen, info);

	return fp;
}
