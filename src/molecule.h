#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include "sachem.h"

#define Q_ATOM_NUMBER           ((int8_t) -'Q')
#define M_ATOM_NUMBER           ((int8_t) -'M')
#define X_ATOM_NUMBER           ((int8_t) -'X')
#define H_ATOM_NUMBER           1
#define C_ATOM_NUMBER           6
#define MAX_ATOM_IDX            INT16_MAX

#define BOND_LIST_BASE_SIZE     16
#define BOND_BLOCK_SIZE         4
#define HBOND_BLOCK_SIZE        2
#define SPECIAL_BLOCK_SIZE      3


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


typedef int16_t AtomIdx;
typedef int16_t BondIdx;
typedef int16_t MolSize;


typedef struct Molecule
{
    int atomCount;
    int bondCount;
    int8_t *restrict atomNumbers;
    uint8_t *restrict atomHydrogens;
    int8_t *restrict atomCharges;
    int8_t *restrict atomMasses;
    uint8_t *restrict atomStereo;
    uint8_t *restrict bondTypes;
    uint8_t *restrict bondStereo;
    bool *restrict restH;

    BondIdx *restrict bondMatrix;
    AtomIdx *restrict bondLists;
    MolSize *restrict bondListSizes;

    AtomIdx (*restrict contains)[2];
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
            bool *extendedRestH = (bool *) palloc0((size_t) atomCount * sizeof(bool));
            memcpy(extendedRestH, restH, (size_t) heavyAtomCount * sizeof(bool));
            restH = extendedRestH;
        }
    }

    int8_t *atomNumbers = (int8_t *) palloc((size_t) atomCount);
    uint8_t *atomHydrogens = (uint8_t *) palloc0((size_t) atomCount);
    uint8_t *bondTypes = (uint8_t *) palloc((size_t) bondCount);
    int8_t *atomCharges = withCharges ? (int8_t *) palloc0((size_t) atomCount) : NULL;
    int8_t *atomMasses = withIsotopes ? (int8_t *) palloc0((size_t) atomCount) : NULL;
    uint8_t *atomStereo = withStereo ? (uint8_t *) palloc0((size_t) atomCount) : NULL;
    uint8_t *bondStereo = withStereo ? (uint8_t *) palloc0((size_t) bondCount) : NULL;

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

    molecule->bondLists = (AtomIdx *) palloc(BOND_LIST_BASE_SIZE * (size_t) atomCount * sizeof(AtomIdx));
    molecule->bondListSizes = (MolSize *) palloc0((size_t) atomCount * sizeof(MolSize));
    molecule->contains = (AtomIdx (*)[2]) palloc((size_t) bondCount * 2 * sizeof(AtomIdx));
    molecule->bondMatrix = (BondIdx *) palloc((size_t) atomCount * (size_t) atomCount * sizeof(BondIdx));


    for(int i = 0; i < xAtomCount; i++)
        atomNumbers[i] = (int8_t) data[i];

    for(int i = xAtomCount; i < heavyAtomCount; i++)
        atomNumbers[i] = C_ATOM_NUMBER;

    for(int i = heavyAtomCount; i < atomCount; i++)
        atomNumbers[i] = H_ATOM_NUMBER;

    data += xAtomCount;


    for(int i = 0; i < atomCount * atomCount; i++)
        molecule->bondMatrix[i] = -1;

    BondIdx boundIdx = 0;

    for(int i = 0; i < xBondCount; i++)
    {
        int offset = i * BOND_BLOCK_SIZE;

        bondTypes[i] = data[offset + 3];

        int b0 = data[offset + 0];
        int b1 = data[offset + 1];
        int b2 = data[offset + 2];

        AtomIdx x = (AtomIdx) (b0 | (b1 << 4 & 0xF00));
        AtomIdx y = (AtomIdx) (b2 | (b1 << 8 & 0xF00));

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


        AtomIdx idx = value & 0xFFF;

        if(idx < atomCount)
            atomHydrogens[idx]++;


        if(extended)
        {
            AtomIdx idy = (AtomIdx) (heavyAtomCount + i);

            if(idx >= heavyAtomCount)
                atomHydrogens[idy]++;

            bondTypes[boundIdx] = (uint8_t ) (data[offset] >> 4);


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
                    atomMasses[idx] = (int8_t) data[offset + 2];
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

    molecule->atomNumbers = (int8_t *) palloc((size_t) atomCount);
    molecule->bondTypes = (uint8_t *) palloc((size_t) bondCount);
    molecule->bondLists = (AtomIdx *) palloc(BOND_LIST_BASE_SIZE * (size_t) atomCount * sizeof(AtomIdx));
    molecule->bondListSizes = (MolSize *) palloc0((size_t) atomCount * sizeof(MolSize));
    molecule->contains = (AtomIdx (*)[2]) palloc((size_t) bondCount * 2 * sizeof(AtomIdx));
    molecule->bondMatrix = (BondIdx *) palloc((size_t) atomCount * (size_t) atomCount * sizeof(BondIdx));


    for(int i = 0; i < xAtomCount; i++)
        molecule->atomNumbers[i] = (int8_t) data[i];

    for(int i = xAtomCount; i < atomCount; i++)
        molecule->atomNumbers[i] = C_ATOM_NUMBER;

    data += xAtomCount;


    for(int i = 0; i < atomCount * atomCount; i++)
        molecule->bondMatrix[i] = -1;

    BondIdx boundIdx = 0;

    for(int i = 0; i < xBondCount; i++)
    {
        int offset = i * BOND_BLOCK_SIZE;

        molecule->bondTypes[i] = data[offset + 3];

        int b0 = data[offset + 0];
        int b1 = data[offset + 1];
        int b2 = data[offset + 2];

        AtomIdx x = (AtomIdx) (b0 | (b1 << 4 & 0xF00));
        AtomIdx y = (AtomIdx) (b2 | (b1 << 8 & 0xF00));

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

        int x = b0 | (b1 << 4 & 0xF00);
        int y = b2 | (b1 << 8 & 0xF00);

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


inline bool molecule_is_pseudo_atom(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule->atomNumbers[atom] < 0;
}


inline int8_t molecule_get_atom_number(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule->atomNumbers[atom];
}


inline bool molecule_is_metal(const Molecule *const restrict molecule, AtomIdx atom)
{
    int8_t number = molecule_get_atom_number(molecule, atom);

    return  (number > 2 && number < 5) || (number > 10 && number < 14) ||
            (number > 18 && number < 32) || (number > 36 && number < 51) ||
            (number > 54 && number < 85) || number > 86;
}


inline bool molecule_is_halogen(const Molecule *const restrict molecule, AtomIdx atom)
{
    int8_t number = molecule_get_atom_number(molecule, atom);

    return number == 9 || number == 17 || number == 35 || number == 53 || number == 85;
}


inline uint8_t molecule_get_hydrogen_count(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule->atomHydrogens[atom];
}


inline int8_t molecule_get_formal_charge(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule->atomCharges[atom];
}


inline int8_t molecule_get_atom_mass(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule->atomMasses[atom];
}


inline AtomIdx *molecule_get_bonded_atom_list(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule->bondLists + atom * BOND_LIST_BASE_SIZE;
}


inline MolSize molecule_get_bonded_atom_list_size(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule->bondListSizes[atom];
}


inline bool molecule_has_restH_flags(const Molecule *const restrict molecule)
{
    return molecule->restH != NULL;
}


inline bool molecule_get_atom_restH_flag(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule->restH[atom];
}


inline AtomIdx *molecule_bond_atoms(const Molecule *const restrict molecule, BondIdx bond)
{
    return molecule->contains[bond];
}


inline bool molecule_bond_contains(const Molecule *const restrict molecule, BondIdx bond, AtomIdx atom)
{
    return molecule->contains[bond][0] == atom || molecule->contains[bond][1] == atom;
}


inline AtomIdx molecule_get_other_bond_atom(const Molecule *const restrict molecule, BondIdx bond, AtomIdx atom)
{
    if(molecule->contains[bond][0] == atom)
        return molecule->contains[bond][1];
    else if(molecule->contains[bond][1] == atom)
        return molecule->contains[bond][0];
    return -1;
}


inline AtomIdx molecule_get_opposite_atom(const Molecule *const restrict molecule, AtomIdx centre, AtomIdx atom)
{
    if(molecule_get_bonded_atom_list_size(molecule, centre) != 2)
        return -1;

    AtomIdx *bonded = molecule_get_bonded_atom_list(molecule, centre);

    if(bonded[0] == atom)
        return bonded[1];
    else if(bonded[1] == atom)
        return bonded[0];
    else
        return -1;
}


inline BondIdx molecule_get_bond(const Molecule *const restrict molecule, AtomIdx i, AtomIdx j)
{
    return molecule->bondMatrix[i * molecule->atomCount + j];
}


inline uint8_t molecule_get_bond_type(const Molecule *const restrict molecule, BondIdx bond)
{
    return molecule->bondTypes[bond];
}


inline uint8_t molecule_get_atom_stereo(const Molecule *const restrict molecule, AtomIdx atom)
{
    return molecule->atomStereo[atom];
}


inline uint8_t molecule_get_bond_stereo(const Molecule *const restrict molecule, BondIdx bond)
{
    return molecule->bondStereo[bond];
}


inline AtomIdx molecule_get_last_stereo_bond_ligand(const Molecule *const restrict molecule, AtomIdx atom, BondIdx bond, AtomIdx ligand)
{
    MolSize listSize = molecule_get_bonded_atom_list_size(molecule, atom);

    if(listSize == 2)
    {
        return MAX_ATOM_IDX;
    }
    else if(listSize == 3)
    {
        AtomIdx other = molecule_get_other_bond_atom(molecule, bond, atom);

        for(int i = 0; i < 3; i++)
        {
            AtomIdx a = molecule_get_bonded_atom_list(molecule, atom)[i];

            if(a != other && a != ligand)
                return a;
        }
    }

    elog(NOTICE, "unexpected branch: %s: %i", __FILE__, __LINE__);
    return MAX_ATOM_IDX;
}


inline AtomIdx molecule_get_last_chiral_ligand(const Molecule *const restrict molecule, AtomIdx centre, AtomIdx *ligands)
{
    MolSize listSize = molecule_get_bonded_atom_list_size(molecule, centre);

    if(listSize == 3)
    {
        return MAX_ATOM_IDX;
    }
    else if(listSize == 4)
    {
        for(int i = 0; i < 4; i++)
        {
            AtomIdx ligand = molecule_get_bonded_atom_list(molecule, centre)[i];
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
            AtomIdx atom = centre;
            AtomIdx bonded = molecule_get_bonded_atom_list(molecule, centre)[i];

            while(true)
            {
                MolSize newListSize = molecule_get_bonded_atom_list_size(molecule, bonded);

                if(newListSize == 3)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        AtomIdx o = molecule_get_bonded_atom_list(molecule, bonded)[j];

                        if(o == atom)
                            continue;

                        if(o != ligands[0] && o != ligands[1] && o != ligands[2])
                            return o;
                    }

                    break;
                }
                else if(newListSize == 2)
                {
                    AtomIdx next = molecule_get_opposite_atom(molecule, bonded, atom);

                    if(molecule_get_bond_type(molecule, molecule_get_bond(molecule, bonded, next)) != BOND_DOUBLE)
                        return MAX_ATOM_IDX;

                    atom = bonded;
                    bonded = next;
                }
                else
                {
                    elog(NOTICE, "unexpected branch: %s: %i", __FILE__, __LINE__);
                    return MAX_ATOM_IDX;
                }
            }
        }
    }

    elog(NOTICE, "unexpected branch: %s: %i", __FILE__, __LINE__);
    return MAX_ATOM_IDX;
}


inline uint8_t normalize_atom_stereo(AtomIdx indexes[4], uint8_t stereo)
{
    uint16_t order = 0;

    for(int i = 0; i < 4; i++)
    {
        uint16_t value = 0;

        for(int j = 0; j < 4; j++)
            if(indexes[j] <= indexes[i])
                value++;

        order = (uint16_t) ((order << 4) + value);
    }

    bool reverse = true;

    const uint16_t validReorder[] = {0x1234, 0x1423, 0x1342, 0x2314, 0x2431, 0x2143, 0x3124, 0x3412, 0x3241, 0x4213, 0x4321, 0x4132};

    for(int i = 0; i < 12; i++)
        if(validReorder[i] == order)
            reverse = false;

    if(reverse)
        return (uint8_t) ~stereo;

    return stereo;
}


inline uint8_t normalize_bond_stereo(AtomIdx indexes[4], uint8_t conformation)
{
    bool reverse = false;

    if(indexes[0] > indexes[1])
        reverse = !reverse;

    if(indexes[2] > indexes[3])
        reverse = !reverse;

    if(reverse)
        return (uint8_t) ~conformation;

    return conformation;
}

#endif /* MOLECULE_H_ */
