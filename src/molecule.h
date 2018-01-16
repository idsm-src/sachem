#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <stdbool.h>
#include <stdint.h>
#include "sachem.h"

#define Q_ATOM_NUMBER   (-'Q')
#define H_ATOM_NUMBER        1
#define C_ATOM_NUMBER        6

#define BOND_LIST_BASE_SIZE 16
#define BOND_BLOCK_SIZE      4
#define HBOND_BLOCK_SIZE     2
#define SPECIAL_BLOCK_SIZE   3


enum BondType
{
    BOND_NONE = 0,
    BOND_SINGLE = 1,
    BOND_DOUBLE = 2,
    BOND_TRIPLE = 3,
    BOND_QUADRUPLE = 4,
    BOND_QUINTUPLE = 5,
    BOND_SEXTUPLE = 6,
    BOND_AROMATIC = 11,
    BOND_SINGLE_OR_DOUBLE = 12,
    BOND_SINGLE_OR_AROMATIC = 13,
    BOND_DOUBLE_OR_AROMATIC = 14,
    BOND_ANY = 15
};


enum SpecialRecordType
{
    RECORD_CHARGE = 0,
    RECORD_ISOTOPE = 1,
    RECORD_TETRAHEDRAL_STEREO = 2,
    RECORD_BOND_STEREO = 3
};


enum TetrahedralStereoType
{
    TETRAHEDRAL_STEREO_NONE = 0,
    TETRAHEDRAL_STEREO_CLOCKWISE = 1,
    TETRAHEDRAL_STEREO_ANTI_CLOCKWISE = 2,
    TETRAHEDRAL_STEREO_UNDEFINED = 3
};


enum BondStereoType
{
    BOND_STEREO_NONE = 0,
    BOND_STEREO_OPPOSITE = 1,
    BOND_STEREO_TOGETHER = 2,
    BOND_STEREO_UNDEFINED = 3
};


typedef struct
{
    int atomCount;
    int bondCount;
    uint8_t *restrict atomNumbers;
    uint8_t *restrict atomHydrogens;
    int8_t *restrict atomCharges;
    uint8_t *restrict atomMasses;
    uint8_t *restrict atomStereo;
    uint8_t *restrict bondTypes;
    uint8_t *restrict bondStereo;
    bool *restrict restH;

    int *restrict bondMatrix;
    int *restrict bondLists;
    int *restrict bondListSizes;

    int (*restrict contains)[2];
} Molecule;


inline void molecule_init(Molecule *const molecule, uint8_t *data, bool *restH)
{
    int xAtomCount = *data << 8 | *(data + 1);
    data += 2;

    int cAtomCount = *data << 8 | *(data + 1);
    data += 2;

    int hAtomCount = *data << 8 | *(data + 1);
    data += 2;

    int xBondCount = *data << 8 | *(data + 1);
    data += 2;

    int specialCount = *data << 8 | *(data + 1);
    data += 2;


    int atomCount = xAtomCount + cAtomCount;
    int bondCount = xBondCount;

    uint8_t *atomNumbers = (uint8_t *) palloc(atomCount);
    uint8_t *atomHydrogens = (uint8_t *) palloc0(atomCount);
    int8_t  *atomCharges = (int8_t *) palloc0(atomCount);
    uint8_t *atomMasses = (uint8_t *) palloc0(atomCount);
    uint8_t *atomStereo = (uint8_t *) palloc0(atomCount);
    uint8_t *bondTypes = (uint8_t *) palloc(bondCount);
    uint8_t *bondStereo = (uint8_t *) palloc0(bondCount);

    molecule->atomCount = atomCount;
    molecule->bondCount = bondCount;
    molecule->atomNumbers = atomNumbers;
    molecule->atomHydrogens = atomHydrogens;
    molecule->atomCharges = atomCharges;
    molecule->atomMasses = atomMasses;
    molecule->atomStereo = atomStereo;
    molecule->bondTypes = bondTypes;
    molecule->bondStereo = bondStereo;
    molecule->restH = restH;

    molecule->bondLists = (int *) palloc(BOND_LIST_BASE_SIZE * atomCount * sizeof(int));
    molecule->bondListSizes = (int *) palloc0(atomCount * sizeof(int));
    molecule->contains = (int (*)[2]) palloc(bondCount * 2 * sizeof(int));
    molecule->bondMatrix = (int *) palloc(atomCount * atomCount * sizeof(int));


    for(int i = 0; i < xAtomCount; i++)
        atomNumbers[i] = data[i];

    for(int i = xAtomCount; i < xAtomCount + cAtomCount; i++)
        atomNumbers[i] = C_ATOM_NUMBER;

    data += xAtomCount;


    for(int i = 0; i < xBondCount; i++)
        bondTypes[i] = data[BOND_BLOCK_SIZE * i + 3];


    for(int i = 0; i < atomCount * atomCount; i++)
        molecule->bondMatrix[i] = -1;

    for(int i = 0; i < xBondCount; i++)
    {
        int offset = i * BOND_BLOCK_SIZE;

        bondTypes[i] = data[offset + 3];

        int b0 = data[offset + 0];
        int b1 = data[offset + 1];
        int b2 = data[offset + 2];

        int x = b0 | b1 << 4 & 0xF00;
        int y = b2 | b1 << 8 & 0xF00;

        if(x >= atomCount || y >= atomCount)
        {
            if(x < atomCount)
                atomHydrogens[x]++;

            if(y < atomCount)
                atomHydrogens[y]++;

            molecule->bondCount--;
            continue;
        }

        molecule->bondLists[x * BOND_LIST_BASE_SIZE + molecule->bondListSizes[x]++] = y;
        molecule->bondLists[y * BOND_LIST_BASE_SIZE + molecule->bondListSizes[y]++] = x;

        if(unlikely(molecule->bondListSizes[x] == BOND_LIST_BASE_SIZE || molecule->bondListSizes[y] == BOND_LIST_BASE_SIZE))
            elog(ERROR, "%s: too high atom valence", __func__);

        molecule->bondMatrix[x * molecule->atomCount + y] = i;
        molecule->bondMatrix[y * molecule->atomCount + x] = i;

        molecule->contains[i][0] = x;
        molecule->contains[i][1] = y;
    }

    data += xBondCount * BOND_BLOCK_SIZE;


    for(int i = 0; i < hAtomCount; i++)
    {
        int offset = i * HBOND_BLOCK_SIZE;
        int value = data[offset + 0] * 256 | data[offset + 1];

        if(value != 0)
        {
            int idx = value & 0xFFF;

            if(idx < atomCount)
                atomHydrogens[idx]++;
        }
    }

    data += hAtomCount * HBOND_BLOCK_SIZE;


    for(int i = 0; i < specialCount; i++)
    {
        int offset = i * SPECIAL_BLOCK_SIZE;
        int value = data[offset + 0] * 256 | data[offset + 1];
        int idx = value & 0xFFF;

        switch(data[offset] >> 4)
        {
            case RECORD_CHARGE:
                if(idx < atomCount)
                    atomCharges[idx] = (int8_t) data[offset + 2];

            case RECORD_ISOTOPE:
                if(idx < atomCount)
                    atomMasses[idx] = data[offset + 2];

            case RECORD_TETRAHEDRAL_STEREO:
                if(idx < atomCount)
                    atomStereo[idx] = data[offset + 2];

            case RECORD_BOND_STEREO:
                if(idx < bondCount)
                    bondStereo[idx] = data[offset + 2];
        }
    }
}


inline bool molecule_is_pseudo_atom(const Molecule *const restrict molecule, int i)
{
    return ((int8_t) molecule->atomNumbers[i]) < 0;
}


inline int8_t molecule_get_atom_number(const Molecule *const restrict molecule, int i)
{
    return molecule->atomNumbers[i];
}


inline uint8_t molecule_get_hydrogen_count(const Molecule *const restrict molecule, int i)
{
    if(!molecule->atomHydrogens)
        return 0;

    return molecule->atomHydrogens[i];
}


inline int8_t molecule_get_formal_charge(const Molecule *const restrict molecule, int i)
{
    if(!molecule->atomCharges)
        return 0;

    return molecule->atomCharges[i];
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


inline uint8_t molecule_get_bond_type(const Molecule *const restrict molecule, int b)
{
    return molecule->bondTypes[b];
}


inline uint8_t molecule_get_atom_stereo(const Molecule *const restrict molecule, int i)
{
    if(likely(!molecule->atomStereo))
        return TETRAHEDRAL_STEREO_NONE;

    return molecule->atomStereo[i];
}


inline uint8_t molecule_get_bond_stereo(const Molecule *const restrict molecule, int i)
{
    if(likely(!molecule->bondStereo))
        return BOND_STEREO_NONE;

    return molecule->bondStereo[i];
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
