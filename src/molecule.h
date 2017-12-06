#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <stdbool.h>
#include <stdint.h>
#include "sachem.h"

#define Q_ATOM_NUMBER   (-'Q')
#define H_ATOM_NUMBER        1
#define C_ATOM_NUMBER        6

#define BOND_LIST_BASE_SIZE 16
#define ATOM_BLOCK_SIZE      3
#define BOND_BLOCK_SIZE      4

#define STEREO_MASK       0x03
#define NONE_STEREO       0x00
#define UNDEF_STEREO      0x03


typedef struct
{
    int atomCount;
    int bondCount;
    uint8_t *restrict atoms;
    uint8_t *restrict bonds;
    bool *restrict restH;

    int *restrict bondMatrix;
    int *restrict bondLists;
    int *restrict bondListSizes;

    int (*restrict contains)[2];
} Molecule;


#ifndef MOLECULE_H_NO_POSTGRES
inline void molecule_init(Molecule *const molecule, int atomsLength, uint8_t *atoms, int bondsLength, uint8_t *bonds, bool *restH, bool explitHonly)
{
    const int atomCount = atomsLength / ATOM_BLOCK_SIZE;
    const int bondCount = bondsLength / BOND_BLOCK_SIZE;

    molecule->atomCount = atomCount;
    molecule->bondCount = bondCount;
    molecule->atoms = atoms;
    molecule->bonds = bonds;
    molecule->restH = restH;

    molecule->bondLists = (int *) palloc(BOND_LIST_BASE_SIZE * atomCount * sizeof(int));
    molecule->bondListSizes = (int *) palloc0(atomCount * sizeof(int));
    molecule->contains = (int (*)[2]) palloc(bondCount * 2 * sizeof(int));
    molecule->bondMatrix = (int *) palloc(atomCount * atomCount * sizeof(int));


    for(int i = 0; i < atomCount; i++)
    {
        int offset = ATOM_BLOCK_SIZE * i + 1;

        if(explitHonly)
            molecule->atoms[offset] = molecule->atoms[offset] >> 4;
        else
            molecule->atoms[offset] = (molecule->atoms[offset] & 0x0F) + (molecule->atoms[offset] >> 4);
    }


    for(int i = 0; i < atomCount * atomCount; i++)
        molecule->bondMatrix[i] = -1;


    for(int i = 0; i < molecule->bondCount; i++)
    {
        int offset = i * BOND_BLOCK_SIZE;

        int b0 = bonds[offset + 0];
        int b1 = bonds[offset + 1];
        int b2 = bonds[offset + 2];

        if(bonds[offset + 3] & 128)
            bonds[offset + 3] = 128;

        int x = b0 | b1 << 4 & 0xF00;
        int y = b2 | b1 << 8 & 0xF00;

        molecule->bondLists[x * BOND_LIST_BASE_SIZE + molecule->bondListSizes[x]++] = y;
        molecule->bondLists[y * BOND_LIST_BASE_SIZE + molecule->bondListSizes[y]++] = x;

        if(unlikely(molecule->bondListSizes[x] == BOND_LIST_BASE_SIZE || molecule->bondListSizes[y] == BOND_LIST_BASE_SIZE))
            elog(ERROR, "%s: too high atom valence", __func__);

        molecule->bondMatrix[x * molecule->atomCount + y] = i;
        molecule->bondMatrix[y * molecule->atomCount + x] = i;

        molecule->contains[i][0] = x;
        molecule->contains[i][1] = y;
    }
}
#endif


inline bool molecule_is_pseudo_atom(const Molecule *const restrict molecule, int i)
{
    int offset = ATOM_BLOCK_SIZE * i;

    return ((int8_t) molecule->atoms[offset]) < 0;
}


inline int8_t molecule_get_atom_number(const Molecule *const restrict molecule, int i)
{
    int offset = ATOM_BLOCK_SIZE * i;

    return molecule->atoms[offset];
}


inline uint8_t molecule_get_hydrogen_count(const Molecule *const restrict molecule, int i)
{
    int offset = ATOM_BLOCK_SIZE * i;

    return molecule->atoms[offset + 1];
}


inline int8_t molecule_get_formal_charge(const Molecule *const restrict molecule, int i)
{
    int offset = ATOM_BLOCK_SIZE * i + 2;

    return ((int8_t) molecule->atoms[offset]) / 4;
}


inline int *molecule_get_bond_list(const Molecule *const restrict molecule, int i)
{
    return molecule->bondLists + i * BOND_LIST_BASE_SIZE;
}


inline int molecule_get_bond_list_size(const Molecule *const restrict molecule, int i)
{
    return molecule->bondListSizes[i];
}


inline bool molecule_get_atom_restH_flag(const Molecule *const restrict molecule, int i)
{
    if(likely(!molecule->restH))
        return false;

    return molecule->restH[i];
}


inline int *molecule_bond_atoms(const Molecule *const restrict molecule, int b)
{
    return molecule->contains[b];
}


inline bool molecule_bond_contains(const Molecule *const restrict molecule, int b, int i)
{
    return molecule->contains[b][0] == i || molecule->contains[b][1] == i;
}


inline int molecule_get_bond_connected_atom(const Molecule *const restrict molecule, int b, int i)
{
    if(molecule->contains[b][0] == i)
        return molecule->contains[b][1];
    else if(molecule->contains[b][1] == i)
        return molecule->contains[b][0];
    return -1;
}


inline int molecule_get_bond(const Molecule *const restrict molecule, int i, int j)
{
    return molecule->bondMatrix[i * molecule->atomCount + j];
}


inline uint8_t molecule_get_bond_data(const Molecule *const restrict molecule, int b)
{
    int offset = BOND_BLOCK_SIZE * b + 3;

    return molecule->bonds[offset] & 0xFC;
}


inline uint8_t molecule_get_atom_stereo(const Molecule *const restrict molecule, int i)
{
    int offset = ATOM_BLOCK_SIZE * i + 2;

    return molecule->atoms[offset] & STEREO_MASK;
}


inline uint8_t molecule_get_bond_stereo(const Molecule *const restrict molecule, int b)
{
    int offset = BOND_BLOCK_SIZE * b + 3;

    return molecule->bonds[offset] & STEREO_MASK;
}


inline uint8_t normalize_atom_stereo(int indexes[4], uint8_t stereo)
{
    uint16_t order = 0;

    for(int i = 0; i < 4; i++)
    {
        int value = 0;

        for(int j = 0; j < 4; j++)
            if(indexes[j] <= indexes[i])
                value++;

        order = (order << 4) + value;
    }

    bool reverse = true;

    const uint16_t validReorder[] = {0x1234, 0x1423, 0x1342, 0x2314, 0x2431, 0x2143, 0x3124, 0x3412, 0x3241, 0x4213, 0x4321, 0x4132};

    for(int i = 0; i < 12; i++)
        if(validReorder[i] == order)
            reverse = false;

    if(reverse)
        return ~stereo;

    return stereo;
}


inline uint8_t normalize_bond_stereo(int indexes[4], uint8_t conformation)
{
    bool reverse = false;

    if(indexes[0] > indexes[1])
        reverse = !reverse;

    if(indexes[2] > indexes[3])
        reverse = !reverse;

    if(reverse)
        return ~conformation;

    return conformation;
}

#endif /* MOLECULE_H_ */
