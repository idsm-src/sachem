#include <list>
#include <vector>
#include <algorithm>
#include "AtomFingerprint.hpp"

extern "C"
{
#include "molecule.h"
}


std::map<uint32_t, int> atom_fingerprint_get(const Molecule *molecule, BitInfo *info)
{
	std::map<uint32_t, int> fp;

	for(int i = 0; i < molecule->atomCount; i++)
	{
	    if(molecule_get_atom_number(molecule, i) <= H_ATOM_NUMBER)
	        continue;

	    uint32_t hsh = molecule_get_atom_number(molecule, i);

        fp[hsh] += 1;

        if(info)
            (*info)[hsh].insert(i);
	}

	return fp;
}
