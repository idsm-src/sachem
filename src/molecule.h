#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include "sachem.h"

#define Q_ATOM_NUMBER   (-'Q')
#define M_ATOM_NUMBER   (-'M')
#define X_ATOM_NUMBER   (-'X')
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


typedef struct Molecule
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


inline void molecule_init(Molecule *const molecule, const uint8_t *data, bool *restH, bool extended,
        bool withCharges, bool withIsotopes, bool withStereo)
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


    int heavyAtomCount = xAtomCount + cAtomCount;
    int atomCount = xAtomCount + cAtomCount;
    int bondCount = xBondCount;

    if(extended)
    {
        atomCount += hAtomCount;
        bondCount += hAtomCount;

        if(restH != NULL)
        {
            bool *extendedRestH = (bool *) palloc0(atomCount * sizeof(bool));
            memcpy(extendedRestH, restH, heavyAtomCount * sizeof(bool));
            restH = extendedRestH;
        }
    }

    uint8_t *atomNumbers = (uint8_t *) palloc(atomCount);
    uint8_t *atomHydrogens = (uint8_t *) palloc0(atomCount);
    uint8_t *bondTypes = (uint8_t *) palloc(bondCount);
    int8_t  *atomCharges = withCharges ? (int8_t *) palloc0(atomCount) : NULL;
    uint8_t *atomMasses = withIsotopes ? (uint8_t *) palloc0(atomCount) : NULL;
    uint8_t *atomStereo = withStereo ? (uint8_t *) palloc0(atomCount) : NULL;
    uint8_t *bondStereo = withStereo ? (uint8_t *) palloc0(bondCount) : NULL;

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

    for(int i = xAtomCount; i < heavyAtomCount; i++)
        atomNumbers[i] = C_ATOM_NUMBER;

    for(int i = heavyAtomCount; i < atomCount; i++)
        atomNumbers[i] = H_ATOM_NUMBER;

    data += xAtomCount;


    for(int i = 0; i < atomCount * atomCount; i++)
        molecule->bondMatrix[i] = -1;

    int boundIdx = 0;

    for(int i = 0; i < xBondCount; i++)
    {
        int offset = i * BOND_BLOCK_SIZE;

        bondTypes[i] = data[offset + 3];

        int b0 = data[offset + 0];
        int b1 = data[offset + 1];
        int b2 = data[offset + 2];

        int x = b0 | b1 << 4 & 0xF00;
        int y = b2 | b1 << 8 & 0xF00;

        if(x >= heavyAtomCount && y < atomCount)
            atomHydrogens[y]++;

        if(y >= heavyAtomCount && x < atomCount)
            atomHydrogens[x]++;

        if(x >= atomCount || y >= atomCount)
        {
            molecule->bondCount--;
            continue;
        }

        molecule->bondLists[x * BOND_LIST_BASE_SIZE + molecule->bondListSizes[x]++] = y;
        molecule->bondLists[y * BOND_LIST_BASE_SIZE + molecule->bondListSizes[y]++] = x;

        if(unlikely(molecule->bondListSizes[x] == BOND_LIST_BASE_SIZE || molecule->bondListSizes[y] == BOND_LIST_BASE_SIZE))
            elog(ERROR, "%s: too high atom valence", __func__);

        molecule->bondMatrix[x * molecule->atomCount + y] = boundIdx;
        molecule->bondMatrix[y * molecule->atomCount + x] = boundIdx;

        molecule->contains[boundIdx][0] = x;
        molecule->contains[boundIdx][1] = y;

        boundIdx++;
    }

    data += xBondCount * BOND_BLOCK_SIZE;


    for(int i = 0; i < hAtomCount; i++)
    {
        int offset = i * HBOND_BLOCK_SIZE;
        int value = data[offset + 0] * 256 | data[offset + 1];

        if(value == 0)
        {
            if(extended)
                molecule->bondCount--;

            continue;
        }


        int idx = value & 0xFFF;

        if(idx < atomCount)
            atomHydrogens[idx]++;


        if(extended)
        {
            int idy = heavyAtomCount + i;

            if(idx >= heavyAtomCount)
                atomHydrogens[idy]++;

            bondTypes[boundIdx] = data[offset] >> 4;


            molecule->bondLists[idx * BOND_LIST_BASE_SIZE + molecule->bondListSizes[idx]++] = idy;
            molecule->bondLists[idy * BOND_LIST_BASE_SIZE + molecule->bondListSizes[idy]++] = idx;

            if(unlikely(molecule->bondListSizes[idx] == BOND_LIST_BASE_SIZE || molecule->bondListSizes[idy] == BOND_LIST_BASE_SIZE))
                elog(ERROR, "%s: too high atom valence", __func__);

            molecule->bondMatrix[idx * molecule->atomCount + idy] = boundIdx;
            molecule->bondMatrix[idy * molecule->atomCount + idx] = boundIdx;

            molecule->contains[boundIdx][0] = idx;
            molecule->contains[boundIdx][1] = idy;

            boundIdx++;
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
                if(withCharges && idx < atomCount)
                    atomCharges[idx] = (int8_t) data[offset + 2];
                break;

            case RECORD_ISOTOPE:
                if(withIsotopes && idx < atomCount)
                    atomMasses[idx] = data[offset + 2];
                break;

            case RECORD_TETRAHEDRAL_STEREO:
                if(withStereo && idx < atomCount)
                    atomStereo[idx] = data[offset + 2];
                break;

            case RECORD_BOND_STEREO:
                if(withStereo && idx < xBondCount) // not for H bond
                    bondStereo[idx] = data[offset + 2];
                break;
        }
    }
}


inline void molecule_simple_free(Molecule *const molecule)
{
    if(molecule->atomNumbers)
        pfree(molecule->atomNumbers);

    if(molecule->bondTypes)
        pfree(molecule->bondTypes);

    if(molecule->bondLists)
        pfree(molecule->bondLists);

    if(molecule->bondListSizes)
        pfree(molecule->bondListSizes);

    if(molecule->contains)
        pfree(molecule->contains);

    if(molecule->bondMatrix)
        pfree(molecule->bondMatrix);
}


inline void molecule_simple_init(Molecule *const molecule, const uint8_t *data)
{
    memset(molecule, 0, sizeof(Molecule));


    int xAtomCount = *data << 8 | *(data + 1);
    data += 2;

    int cAtomCount = *data << 8 | *(data + 1);
    data += 4;

    int xBondCount = *data << 8 | *(data + 1);
    data += 4;


    int atomCount = xAtomCount + cAtomCount;
    int bondCount = xBondCount;


    molecule->atomCount = atomCount;
    molecule->bondCount = bondCount;

    molecule->atomNumbers = (uint8_t *) palloc(atomCount);
    molecule->bondTypes = (uint8_t *) palloc(bondCount);
    molecule->bondLists = (int *) palloc(BOND_LIST_BASE_SIZE * atomCount * sizeof(int));
    molecule->bondListSizes = (int *) palloc0(atomCount * sizeof(int));
    molecule->contains = (int (*)[2]) palloc(bondCount * 2 * sizeof(int));
    molecule->bondMatrix = (int *) palloc(atomCount * atomCount * sizeof(int));


    for(int i = 0; i < xAtomCount; i++)
        molecule->atomNumbers[i] = data[i];

    for(int i = xAtomCount; i < atomCount; i++)
        molecule->atomNumbers[i] = C_ATOM_NUMBER;

    data += xAtomCount;


    for(int i = 0; i < atomCount * atomCount; i++)
        molecule->bondMatrix[i] = -1;

    int boundIdx = 0;

    for(int i = 0; i < xBondCount; i++)
    {
        int offset = i * BOND_BLOCK_SIZE;

        molecule->bondTypes[i] = data[offset + 3];

        int b0 = data[offset + 0];
        int b1 = data[offset + 1];
        int b2 = data[offset + 2];

        int x = b0 | b1 << 4 & 0xF00;
        int y = b2 | b1 << 8 & 0xF00;

        if(x >= atomCount || y >= atomCount)
        {
            molecule->bondCount--;
            continue;
        }

        molecule->bondLists[x * BOND_LIST_BASE_SIZE + molecule->bondListSizes[x]++] = y;
        molecule->bondLists[y * BOND_LIST_BASE_SIZE + molecule->bondListSizes[y]++] = x;

        if(unlikely(molecule->bondListSizes[x] == BOND_LIST_BASE_SIZE || molecule->bondListSizes[y] == BOND_LIST_BASE_SIZE))
            elog(ERROR, "%s: too high atom valence", __func__);

        molecule->bondMatrix[x * molecule->atomCount + y] = boundIdx;
        molecule->bondMatrix[y * molecule->atomCount + x] = boundIdx;

        molecule->contains[boundIdx][0] = x;
        molecule->contains[boundIdx][1] = y;

        boundIdx++;
    }
}


inline bool molecule_is_extended_search_needed(uint8_t *data, bool withCharges, bool withIsotopes)
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

    int heavyAtomCount = xAtomCount + cAtomCount;


    for(int i = 0; i < xAtomCount; i++)
        if(((int8_t) data[i]) < 0)
            return true;

    data += xAtomCount;


    int hBonds[hAtomCount];

    for(int i = 0; i < hAtomCount; i++)
        hBonds[i] = 0;

    for(int i = 0; i < xBondCount; i++)
    {
        int offset = i * BOND_BLOCK_SIZE;

        int b0 = data[offset + 0];
        int b1 = data[offset + 1];
        int b2 = data[offset + 2];

        int x = b0 | b1 << 4 & 0xF00;
        int y = b2 | b1 << 8 & 0xF00;

        if(x >= heavyAtomCount)
            hBonds[x - heavyAtomCount]++;

        if(y >= heavyAtomCount)
            hBonds[y - heavyAtomCount]++;
    }

    data += xBondCount * BOND_BLOCK_SIZE;


    for(int i = 0; i < hAtomCount; i++)
    {
        int offset = i * HBOND_BLOCK_SIZE;
        int value = data[offset + 0] * 256 | data[offset + 1];

        if(value == 0)
            continue;

        hBonds[i]++;

        int idx = value & 0xFFF;

        if(idx >= heavyAtomCount)
            hBonds[idx - heavyAtomCount]++;
    }

    data += hAtomCount * HBOND_BLOCK_SIZE;


    for(int i = 0; i < hAtomCount; i++)
        if(hBonds[i] != 1)
            return true;


    if(!withCharges && !withIsotopes)
        return false;

    for(int i = 0; i < specialCount; i++)
    {
        int offset = i * SPECIAL_BLOCK_SIZE;
        int value = data[offset + 0] * 256 | data[offset + 1];
        int idx = value & 0xFFF;

        switch(data[offset] >> 4)
        {
            case RECORD_CHARGE:
                if(withCharges && idx >= heavyAtomCount)
                    return true;
                break;

            case RECORD_ISOTOPE:
                if(withIsotopes && idx >= heavyAtomCount)
                    return true;
                break;
        }
    }

    return false;
}


inline bool molecule_is_pseudo_atom(const Molecule *const restrict molecule, int i)
{
    return ((int8_t) molecule->atomNumbers[i]) < 0;
}


inline int8_t molecule_get_atom_number(const Molecule *const restrict molecule, int i)
{
    return molecule->atomNumbers[i];
}


inline bool molecule_is_metal(const Molecule *const restrict molecule, int i)
{
    int number = molecule_get_atom_number(molecule, i);

    return  (number > 2 && number < 5) || (number > 10 && number < 14) ||
            (number > 18 && number < 32) || (number > 36 && number < 51) ||
            (number > 54 && number < 85) || number > 86;
}


inline bool molecule_is_halogen(const Molecule *const restrict molecule, int i)
{
    int number = molecule_get_atom_number(molecule, i);

    return number == 9 || number == 17 || number == 35 || number == 53 || number == 85;
}


inline uint8_t molecule_get_hydrogen_count(const Molecule *const restrict molecule, int i)
{
    return molecule->atomHydrogens[i];
}


inline int8_t molecule_get_formal_charge(const Molecule *const restrict molecule, int i)
{
    return molecule->atomCharges[i];
}


inline int8_t molecule_get_atom_mass(const Molecule *const restrict molecule, int i)
{
    return molecule->atomMasses[i];
}


inline int *molecule_get_bond_list(const Molecule *const restrict molecule, int i)
{
    return molecule->bondLists + i * BOND_LIST_BASE_SIZE;
}


inline int molecule_get_bond_list_size(const Molecule *const restrict molecule, int i)
{
    return molecule->bondListSizes[i];
}


inline bool molecule_has_restH_flags(const Molecule *const restrict molecule)
{
    return molecule->restH != NULL;
}


inline bool molecule_get_atom_restH_flag(const Molecule *const restrict molecule, int i)
{
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


inline int molecule_get_atom_connected_bond(const Molecule *const restrict molecule, int i, int b)
{
    if(molecule_get_bond_list_size(molecule, i) != 2)
        return -1;

    int *bonds = molecule_get_bond_list(molecule, i);

    if(bonds[0] == b)
        return bonds[1];
    else if(bonds[1] == b)
        return bonds[0];
    else
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
    return molecule->atomStereo[i];
}


inline uint8_t molecule_get_bond_stereo(const Molecule *const restrict molecule, int i)
{
    return molecule->bondStereo[i];
}


inline int molecule_get_last_stereo_bond_ligand(const Molecule *const restrict molecule, int atom, int bond, int ligand)
{
    int listSize = molecule_get_bond_list_size(molecule, atom);

    if(listSize == 2)
    {
        return INT_MAX;
    }
    else if(listSize == 3)
    {
        for(int i = 0; i < 3; i++)
        {
            int b = molecule_get_bond_list(molecule, atom)[i];

            if(b == bond)
                continue;

            int a = molecule_get_bond_connected_atom(molecule, b, atom);

            if(a == ligand)
                continue;

            return a;
        }
    }

    elog(NOTICE, "unexpected branch: %s: %i", __FILE__, __LINE__);
    return INT_MAX;
}


inline int molecule_get_last_chiral_ligand(const Molecule *const restrict molecule, int centre, int *ligands)
{
    int listSize = molecule_get_bond_list_size(molecule, centre);

    if(listSize == 3)
    {
        return INT_MAX;
    }
    else if(listSize == 4)
    {
        for(int i = 0; i < 4; i++)
        {
            int ligand = molecule_get_bond_list(molecule, centre)[i];
            bool contains = false;

            for(int j = 0; j < 3; j++)
                if(ligands[j] == ligand)
                    contains = true;

            if(!contains)
                return ligand;
        }
    }
    else if(listSize == 2)
    {
        for(int i = 0; i < 2; i++)
        {
            int atom = centre;
            int bond = molecule_get_bond_list(molecule, centre)[i];

            while(true)
            {
                atom = molecule_get_bond_connected_atom(molecule, bond, atom);
                int newListSize = molecule_get_bond_list_size(molecule, atom);

                if(newListSize == 3)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        int b = molecule_get_bond_list(molecule, atom)[j];

                        if(b == bond)
                            continue;

                        int a = molecule_get_bond_connected_atom(molecule, b, atom);

                        if(a != ligands[0] && a != ligands[1] && a != ligands[2])
                            return a;
                    }

                    break;
                }
                else if(newListSize == 2)
                {
                    bond = molecule_get_atom_connected_bond(molecule, atom, bond);

                    if(molecule_get_bond_type(molecule, bond) != BOND_DOUBLE)
                        return INT_MAX;
                }
                else
                {
                    elog(NOTICE, "unexpected branch: %s: %i", __FILE__, __LINE__);
                    return INT_MAX;
                }
            }
        }
    }

    elog(NOTICE, "unexpected branch: %s: %i", __FILE__, __LINE__);
    return INT_MAX;
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
