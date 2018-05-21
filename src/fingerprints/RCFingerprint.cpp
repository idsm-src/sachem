#include <set>
#include <list>
#include <vector>
#include <algorithm>
#include "RCFingerprint.hpp"

extern "C"
{
#include "molecule.h"
}


struct AtomDesc
{
    uint32_t hash;
    AtomIdx atom;
    std::set<BondIdx> cover;

    AtomDesc() {}

    AtomDesc(uint32_t hash, AtomIdx atom) : hash(hash), atom(atom) {}

    AtomDesc(uint32_t hash, AtomIdx atom, std::set<BondIdx> &&cover) : hash(hash), atom(atom), cover(std::move(cover)) {}
};


static inline void update_seed(uint32_t x, uint32_t &seed)
{
    seed ^= (uint32_t) x * 2654435761 + 2654435769 + (seed << 6) + (seed >> 2);
}


static inline uint32_t hash_2(uint32_t a, uint32_t b)
{
    uint32_t seed = 0;
    update_seed(a, seed);
    update_seed(b, seed);
    return seed;
}


static inline uint32_t hash_3(uint32_t a, uint32_t b, uint32_t c)
{
    uint32_t seed = 0;
    update_seed(a, seed);
    update_seed(b, seed);
    update_seed(c, seed);
    return seed;
}


static inline uint32_t hashlist_2_sort(uint32_t a, uint32_t b, std::vector<uint32_t> &l)
{
    uint32_t seed = 0;
    update_seed(a, seed);
    update_seed(b, seed);

    std::sort(l.begin(), l.end());

    for(auto&i : l)
        update_seed(i, seed);

    return seed;
}


static inline bool is_dummy_bond(const Molecule *const restrict molecule, BondIdx bond)
{
    uint8_t type = molecule_get_bond_type(molecule, bond);

    return type == BOND_NONE || type > BOND_AROMATIC;
}


static inline uint32_t bond_hash(const Molecule *const restrict molecule, BondIdx bond)
{
    return molecule_get_bond_type(molecule, bond);
}


static inline bool is_dummy_atom(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule_is_pseudo_atom(molecule, atom);
}


static inline uint32_t atom_hash(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule_get_atom_number(molecule, atom);
}


std::map<uint32_t, int> rc_fingerprint_get(const Molecule *molecule, int minRadius, int maxRadius, BitInfo *info)
{
    std::map<AtomIdx, AtomDesc> desc;
    std::list<AtomDesc> result;


    for(AtomIdx a = 0; a < molecule->atomCount; a++)
    {
        if(molecule_get_atom_number(molecule, a) == H_ATOM_NUMBER || is_dummy_atom(molecule, a))
            continue;

        desc[a] = AtomDesc(atom_hash(molecule, a), a);

        if(minRadius == 0)
            result.push_back(desc[a]);
    }


    for(int radius = 1; radius <= maxRadius; radius++)
    {
        std::map<AtomIdx, AtomDesc> newdesc;
        int updates = 0;

        for(AtomIdx i = 0; i < molecule->atomCount; i++)
        {
            if(molecule_get_atom_number(molecule, i) == H_ATOM_NUMBER || !desc.count(i))
                continue;

            std::vector<uint32_t> hs;
            hs.reserve(molecule->atomCount);
            std::set<BondIdx> newcover;
            bool acceptable = true;

            AtomIdx *neighbors = molecule_get_bonded_atom_list(molecule, i);
            MolSize length = molecule_get_bonded_atom_list_size(molecule, i);

            for(MolSize k = 0; k < length; k++)
            {
                AtomIdx a = neighbors[k];
                BondIdx b = molecule_get_bond(molecule, i, a);

                if(molecule_get_atom_number(molecule, a) == H_ATOM_NUMBER)
                    continue;

                if(is_dummy_atom(molecule, a) || is_dummy_bond(molecule, b) || !desc.count(a))
                {
                    acceptable = false;
                    break;
                }

                newcover.insert(b);
                newcover.insert(desc[a].cover.begin(), desc[a].cover.end());
                hs.push_back(hash_2(bond_hash(molecule, b), desc[a].hash));
            }

            if(!acceptable || newcover.empty())
                continue;


            if(newcover == desc[i].cover)
            {
                newdesc[i] = desc[i];
            }
            else
            {
                newdesc[i] = AtomDesc(hashlist_2_sort(desc[i].hash, newcover.size(), hs), i, std::move(newcover));

                if(radius >= minRadius)
                    result.push_back(newdesc[i]);

                updates++;
            }
        }

        desc.swap(newdesc);

        if(updates == 0)
            break;
    }


    std::map<uint32_t, int> fp;

    for(AtomDesc &i : result)
    {
        fp[i.hash] += 1;

        if(!info)
            continue;

        (*info)[i.hash].insert(i.atom);

        for(BondIdx b : i.cover)
        {
            (*info)[i.hash].insert(molecule_bond_atoms(molecule, b)[0]);
            (*info)[i.hash].insert(molecule_bond_atoms(molecule, b)[1]);
        }
    }

    return fp;
}
