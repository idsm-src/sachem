#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include "memory.h"
#include "sachem.h"


#define UNKNOWN_ATOM_NUMBER     ((int8_t) -'?')
#define R_ATOM_NUMBER           ((int8_t) -'R')
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
    RECORD_BOND_STEREO = 3,
    RECORD_RADICAL = 4,
    RECORD_SGROUP_COUNT = 5
};


enum VariableLengthRecordType
{
    VAR_RECORD_SGROUP = 1,
    VAR_RECORD_ATOM_LABEL = 2
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


enum SgroupType
{
    SGROUP_TYPE_NONE = 0,
    SGROUP_TYPE_SRU = 1,
    SGROUP_TYPE_MOD = 2,
    SGROUP_TYPE_MON = 3,
    SGROUP_TYPE_COP = 4,
    SGROUP_TYPE_GEN = 5,
    SGROUP_TYPE_ANY = 6,
    SGROUP_TYPE_CRO = 7,
    SGROUP_TYPE_MER = 8,
    SGROUP_TYPE_GRA = 9,
    SGROUP_TYPE_COM = 10,
    SGROUP_TYPE_FOR = 11,
    SGROUP_TYPE_MIX = 12
};


enum SgroupSubtype
{
    SGROUP_SUBTYPE_NONE = 0,
    SGROUP_SUBTYPE_ALT = 1,
    SGROUP_SUBTYPE_RAN = 2,
    SGROUP_SUBTYPE_BLO = 3,
    SGROUP_SUBTYPE_UNKNOWN = -1
};


enum SgroupConnectivity
{
    SGROUP_CONNECTIVITY_NONE = 0,
    SGROUP_CONNECTIVITY_HH = 1,
    SGROUP_CONNECTIVITY_HT = 2,
    SGROUP_CONNECTIVITY_EU = 3,
    SGROUP_CONNECTIVITY_UNKNOWN = -1
};


typedef int16_t AtomIdx;
typedef int16_t BondIdx;
typedef int16_t MolSize;


typedef struct
{
    int8_t type;
    int8_t subtype;
    int8_t connectivity;

    int atomCount;
    int bondCount;

    AtomIdx *restrict atoms;
    AtomIdx (*restrict bonds)[2];
}
SGroup;


typedef struct
{
    AtomIdx atom;
    int32_t length;
    uint8_t *label;
}
AtomLabel;


typedef struct
{
    int atomCount;
    int bondCount;

    int heavyAtomCount;
    int heavyBondCount;

    int hydrogenAtomCount;
    int hydrogenBondCount;

    int sgroupCount;
    int labelCount;

    bool extended;

    int8_t *restrict atomNumbers;
    uint8_t *restrict atomHydrogens;
    int8_t *restrict atomCharges;
    int8_t *restrict atomMasses;
    int8_t *restrict atomRadicalTypes;
    uint8_t *restrict atomStereo;
    uint8_t *restrict bondTypes;
    uint8_t *restrict bondStereo;
    uint8_t *restrict restH;
    SGroup *restrict sgroups;
    AtomLabel *restrict labels;

    BondIdx *restrict bondMatrix;
    AtomIdx *restrict bondLists;
    MolSize *restrict bondListSizes;

    AtomIdx (*restrict contains)[2];
}
Molecule;


static inline size_t molecule_mem_size(const uint8_t *restrict data, const uint8_t *restrict restData, bool extended,
        bool withCharges, bool withIsotopes, bool withRadicals, bool withStereo, bool withSGroup,
        bool ignoreChargedHydrogens, bool ignoreHydrogenIsotopes, bool ignoreHydrogenRadicals)
{
    int xAtomCount = data[0] << 8 | data[1];
    int cAtomCount = data[2] << 8 | data[3];
    int hAtomCount = data[4] << 8 | data[5];
    int xBondCount = data[6] << 8 | data[7];
    int specialCount = data[8] << 8 | data[9];

    int atomCount = xAtomCount + cAtomCount;
    int bondCount = xBondCount;

    if(extended)
    {
        atomCount += hAtomCount;
        bondCount += hAtomCount;
    }

    size_t memsize = align_size(sizeof(Molecule));

    if(restData != NULL)
        memsize += align_size(atomCount * sizeof(uint8_t));

    memsize += (2 + withCharges + withIsotopes + withRadicals + withStereo) * align_size(atomCount);
    memsize += (1 + withStereo) * align_size(bondCount);

    memsize += align_size(BOND_LIST_BASE_SIZE * atomCount * sizeof(AtomIdx));
    memsize += align_size(atomCount * sizeof(MolSize));
    memsize += align_size(bondCount * sizeof(AtomIdx[2]));
    memsize += align_size(atomCount * atomCount * sizeof(BondIdx));

    if(!extended && (ignoreChargedHydrogens || ignoreHydrogenIsotopes || ignoreHydrogenRadicals))
        memsize += align_size(sizeof(bool) * hAtomCount);


    data += 10;

    int labelCount = 0;
    int sgroupCount = 0;

    for(int i = 0; i < xAtomCount; i++)
        if((int8_t) data[i] == UNKNOWN_ATOM_NUMBER)
            labelCount++;

    data += xAtomCount + xBondCount * BOND_BLOCK_SIZE + hAtomCount * HBOND_BLOCK_SIZE;

    if(withSGroup)
        for(int i = 0; i < specialCount; i++)
            if(data[i * SPECIAL_BLOCK_SIZE] >> 4 == RECORD_SGROUP_COUNT)
                sgroupCount = (data[i * SPECIAL_BLOCK_SIZE + 0] * 256 | data[i * SPECIAL_BLOCK_SIZE + 1]) & 0xFFF;

    data += specialCount * SPECIAL_BLOCK_SIZE;


    memsize += align_size(labelCount * sizeof(AtomLabel));

    for(int i = 0; i < labelCount; i++)
    {
        int size = data[0] << 24 | data[1] << 16 | data[2] << 8 | data[3];
        data += size;
        memsize += align_size(size - 7);
    }


    memsize += align_size(sgroupCount * sizeof(SGroup));

    for(int i = 0; i < sgroupCount; i++)
    {
        int size = data[0] << 24 | data[1] << 16 | data[2] << 8 | data[3];
        int atomLength = data[8] * 256 | data[9];
        int bondLength = data[10] * 256 | data[11];
        data += size;

        memsize += align_size(atomLength * sizeof(AtomIdx));
        memsize += align_size(bondLength * 2 * sizeof(AtomIdx));
    }


    return memsize;
}


static inline size_t molecule_extended_mem_size(const Molecule *restrict molecule)
{
    int atomCount = molecule->heavyAtomCount + molecule->hydrogenAtomCount;
    int bondCount = molecule->heavyBondCount + molecule->hydrogenBondCount;
    bool withCharges = molecule->atomCharges;
    bool withIsotopes = molecule->atomMasses;
    bool withRadicals = molecule->atomRadicalTypes;
    bool withStereo = molecule->atomStereo;

    size_t memsize = align_size(sizeof(Molecule));

    if(molecule->restH != NULL)
        memsize += align_size(atomCount * sizeof(uint8_t));

    memsize += (2 + withCharges + withIsotopes + withRadicals + withStereo) * align_size(atomCount);
    memsize += (1 + withStereo) * align_size(bondCount);

    memsize += align_size(BOND_LIST_BASE_SIZE * atomCount * sizeof(AtomIdx));
    memsize += align_size(atomCount * sizeof(MolSize));
    memsize += align_size(bondCount * 2 * sizeof(AtomIdx));
    memsize += align_size(atomCount * atomCount * sizeof(BondIdx));


    memsize += align_size(molecule->labelCount * sizeof(AtomLabel));

    for(int i = 0; i < molecule->labelCount; i++)
        memsize += align_size(molecule->labels[i].length);


    memsize += align_size(molecule->sgroupCount * sizeof(SGroup));

    for(int i = 0; i < molecule->sgroupCount; i++)
    {
        memsize += align_size(molecule->sgroups[i].atomCount * sizeof(AtomIdx));
        memsize += align_size(molecule->sgroups[i].bondCount * 2 * sizeof(AtomIdx));
    }

    return memsize;
}


static inline Molecule *molecule_create(void *memory, const uint8_t *restrict data, uint8_t *restrict restData,
        bool extended, bool withCharges, bool withIsotopes, bool withRadicals, bool withStereo, bool withSGroup,
        bool ignoreChargedHydrogens, bool ignoreHydrogenIsotopes, bool ignoreHydrogenRadicals)
{
    ignoreChargedHydrogens &= !extended;
    ignoreHydrogenIsotopes &= !extended;
    ignoreHydrogenRadicals &= !extended;


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
    int heavyBondCount = xBondCount;
    int hydrogenAtomCount = hAtomCount;
    int hydrogenBondCount = hAtomCount;
    int atomCount = xAtomCount + cAtomCount;
    int bondCount = xBondCount;
    int sgroupCount = 0;
    int labelCount = 0;


    if(extended)
    {
        atomCount += hAtomCount;
        bondCount += hAtomCount;
    }


    uint8_t *restrict restH = NULL;

    if(restData != NULL)
    {
        restH = (uint8_t *) alloc_memory(&memory, atomCount * sizeof(uint8_t));
        memcpy(restH, restData, atomCount * sizeof(uint8_t));
    }


    int8_t *restrict atomNumbers = (int8_t *) alloc_memory(&memory, atomCount);
    uint8_t *restrict atomHydrogens = (uint8_t *) alloc_memory_zero(&memory, atomCount);
    uint8_t *restrict bondTypes = (uint8_t *) alloc_memory(&memory, bondCount);
    int8_t *restrict atomCharges = withCharges ? (int8_t *) alloc_memory_zero(&memory, atomCount) : NULL;
    int8_t *restrict atomMasses = withIsotopes ? (int8_t *) alloc_memory_zero(&memory, atomCount) : NULL;
    int8_t *restrict atomRadicalTypes = withRadicals ? (int8_t *) alloc_memory_zero(&memory, atomCount) : NULL;
    uint8_t *restrict atomStereo = withStereo ? (uint8_t *) alloc_memory_zero(&memory, atomCount) : NULL;
    uint8_t *restrict bondStereo = withStereo ? (uint8_t *) alloc_memory_zero(&memory, bondCount) : NULL;

    AtomIdx *restrict bondLists = (AtomIdx *) alloc_memory(&memory, BOND_LIST_BASE_SIZE * atomCount * sizeof(AtomIdx));
    MolSize *restrict bondListSizes = (MolSize *) alloc_memory_zero(&memory, atomCount * sizeof(MolSize));
    AtomIdx (*restrict contains)[2] = (AtomIdx (*)[2]) alloc_memory(&memory, bondCount * sizeof(AtomIdx[2]));
    BondIdx *restrict bondMatrix = (BondIdx *) alloc_memory_one(&memory, atomCount * atomCount * sizeof(BondIdx));


    for(int i = 0; i < xAtomCount; i++)
        atomNumbers[i] = (int8_t) data[i];

    for(int i = xAtomCount; i < heavyAtomCount; i++)
        atomNumbers[i] = C_ATOM_NUMBER;

    for(int i = heavyAtomCount; i < atomCount; i++)
        atomNumbers[i] = H_ATOM_NUMBER;

    for(int i = 0; i < xAtomCount; i++)
        if(atomNumbers[i] == UNKNOWN_ATOM_NUMBER)
            labelCount++;

    data += xAtomCount;


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

        if(x >= heavyAtomCount || y >= heavyAtomCount)
        {
            heavyBondCount--;
            hydrogenBondCount++;
        }

        if(x >= atomCount || y >= atomCount)
        {
            bondCount--;
            continue;
        }

        bondLists[x * BOND_LIST_BASE_SIZE + bondListSizes[x]++] = y;
        bondLists[y * BOND_LIST_BASE_SIZE + bondListSizes[y]++] = x;

        if(unlikely(bondListSizes[x] == BOND_LIST_BASE_SIZE || bondListSizes[y] == BOND_LIST_BASE_SIZE))
            return NULL;

        bondMatrix[x * atomCount + y] = boundIdx;
        bondMatrix[y * atomCount + x] = boundIdx;

        contains[boundIdx][0] = x;
        contains[boundIdx][1] = y;

        boundIdx++;
    }

    data += xBondCount * BOND_BLOCK_SIZE;


    bool *restrict ignoredHydrogen = NULL;;

    if(ignoreChargedHydrogens || ignoreHydrogenIsotopes || ignoreHydrogenRadicals)
    {
        const uint8_t *sdata = data + hAtomCount * HBOND_BLOCK_SIZE;
        ignoredHydrogen = (bool *) alloc_memory_zero(&memory, sizeof(bool) * hAtomCount);

        for(int i = 0; i < specialCount; i++)
        {
            int offset = i * SPECIAL_BLOCK_SIZE;
            int value = sdata[offset + 0] * 256 | sdata[offset + 1];
            int idx = value & 0xFFF;

            switch(sdata[offset] >> 4)
            {
                case RECORD_CHARGE:
                    if(ignoreChargedHydrogens && idx >= heavyAtomCount)
                        ignoredHydrogen[idx - heavyAtomCount] = true;
                    break;

                case RECORD_ISOTOPE:
                    if(ignoreHydrogenIsotopes && idx >= heavyAtomCount)
                        ignoredHydrogen[idx - heavyAtomCount] = true;
                    break;

                case RECORD_RADICAL:
                    if(ignoreHydrogenRadicals && idx >= heavyAtomCount)
                        ignoredHydrogen[idx - heavyAtomCount] = true;
                    break;
            }
        }
    }


    for(int i = 0; i < hAtomCount; i++)
    {
        int offset = i * HBOND_BLOCK_SIZE;
        int value = data[offset + 0] * 256 | data[offset + 1];

        if(value == 0)
        {
            hydrogenBondCount--;

            if(extended)
                bondCount--;

            continue;
        }


        AtomIdx idx = value & 0xFFF;

        if(idx < atomCount && (!(ignoreChargedHydrogens || ignoreHydrogenIsotopes || ignoreHydrogenRadicals) || ignoredHydrogen[i] == false))
            atomHydrogens[idx]++;


        if(extended)
        {
            AtomIdx idy = (AtomIdx) (heavyAtomCount + i);

            if(idx >= heavyAtomCount)
                atomHydrogens[idy]++;

            bondTypes[boundIdx] = (uint8_t ) (data[offset] >> 4);


            bondLists[idx * BOND_LIST_BASE_SIZE + bondListSizes[idx]++] = idy;
            bondLists[idy * BOND_LIST_BASE_SIZE + bondListSizes[idy]++] = idx;

            if(unlikely(bondListSizes[idx] == BOND_LIST_BASE_SIZE || bondListSizes[idy] == BOND_LIST_BASE_SIZE))
                return NULL;

            bondMatrix[idx * atomCount + idy] = boundIdx;
            bondMatrix[idy * atomCount + idx] = boundIdx;

            contains[boundIdx][0] = idx;
            contains[boundIdx][1] = idy;

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

            case RECORD_RADICAL:
                if(withRadicals && idx < atomCount)
                    atomRadicalTypes[idx] = (int8_t) data[offset + 2];
                break;

            case RECORD_TETRAHEDRAL_STEREO:
                if(withStereo && idx < atomCount)
                    atomStereo[idx] = data[offset + 2];
                break;

            case RECORD_BOND_STEREO:
                if(withStereo && idx < xBondCount) // not for H bond
                    bondStereo[idx] = data[offset + 2];
                break;

            case RECORD_SGROUP_COUNT:
                if(withSGroup)
                    sgroupCount = idx;
                break;
        }
    }

    data += specialCount * SPECIAL_BLOCK_SIZE;


    AtomLabel *labels = NULL;

    if(labelCount > 0)
        labels = (AtomLabel *) alloc_memory_zero(&memory, labelCount * sizeof(AtomLabel));

    for(int i = 0; i < labelCount; i++)
    {
        int size = data[0] << 24 | data[1] << 16 | data[2] << 8 | data[3];

        labels[i].atom = data[5] * 256 | data[6];
        labels[i].length = size - 7;
        labels[i].label = (uint8_t *) alloc_memory_zero(&memory, size - 7);
        memcpy(labels[i].label, data + 7, size - 7);
        data += size;
    }


    SGroup *sgroups = NULL;

    if(sgroupCount > 0)
    {
        sgroups = (SGroup *) alloc_memory_zero(&memory, sgroupCount * sizeof(SGroup));

        for(int i = 0; i < sgroupCount; i++)
        {
            sgroups[i].type = data[5];
            sgroups[i].subtype = data[6];
            sgroups[i].connectivity = data[7];

            sgroups[i].atomCount = data[8] * 256 | data[9];
            sgroups[i].bondCount = data[10] * 256 | data[11];

            data += 12;

            sgroups[i].atoms = (AtomIdx *) alloc_memory_zero(&memory, sgroups[i].atomCount * sizeof(AtomIdx));
            sgroups[i].bonds = (AtomIdx (*)[2]) alloc_memory_zero(&memory, sgroups[i].bondCount * 2 * sizeof(AtomIdx));

            for(int a = 0; a < sgroups[i].atomCount; a++)
            {
                sgroups[i].atoms[a] = data[0] * 256 | data[1];
                data += 2;
            }

            for(int b = 0; b < sgroups[i].bondCount; b++)
            {
                sgroups[i].bonds[b][0] = data[0] * 256 | data[1];
                sgroups[i].bonds[b][1] = data[2] * 256 | data[3];
                data += 4;
            }
        }
    }


    Molecule *restrict molecule = alloc_memory(&memory, sizeof(Molecule));

    molecule->atomCount = atomCount;
    molecule->bondCount = bondCount;
    molecule->heavyAtomCount = heavyAtomCount;
    molecule->heavyBondCount = heavyBondCount;
    molecule->hydrogenAtomCount = hydrogenAtomCount;
    molecule->hydrogenBondCount = hydrogenBondCount;
    molecule->sgroupCount = sgroupCount;
    molecule->labelCount = labelCount;
    molecule->extended = extended;
    molecule->atomNumbers = atomNumbers;
    molecule->atomHydrogens = atomHydrogens;
    molecule->atomCharges = atomCharges;
    molecule->atomMasses = atomMasses;
    molecule->atomRadicalTypes = atomRadicalTypes;
    molecule->atomStereo = atomStereo;
    molecule->bondTypes = bondTypes;
    molecule->bondStereo = bondStereo;
    molecule->restH = restH;
    molecule->sgroups = sgroups;
    molecule->labels = labels;
    molecule->bondLists = bondLists;
    molecule->bondListSizes = bondListSizes;
    molecule->contains = contains;
    molecule->bondMatrix = bondMatrix;

    return molecule;
}


static inline Molecule *molecule_extend(void *memory, const Molecule *restrict template)
{
    Molecule *restrict molecule = alloc_memory(&memory, sizeof(Molecule));

    molecule->atomCount = template->heavyAtomCount + template->hydrogenAtomCount;
    molecule->bondCount = template->heavyBondCount + template->hydrogenBondCount;
    molecule->heavyAtomCount = template->heavyAtomCount;
    molecule->heavyBondCount = template->heavyBondCount;
    molecule->hydrogenAtomCount = template->hydrogenAtomCount;
    molecule->hydrogenBondCount = template->hydrogenBondCount;
    molecule->extended = true;

    if(template->restH != NULL)
    {
        molecule->restH = (uint8_t *) alloc_memory(&memory, molecule->atomCount);
        memcpy(molecule->restH, template->restH, template->atomCount);
        memset(molecule->restH + template->atomCount, 0, template->hydrogenAtomCount);
    }
    else
    {
        molecule->restH = NULL;
    }

    molecule->atomNumbers = (int8_t *) alloc_memory(&memory, molecule->atomCount);
    memcpy(molecule->atomNumbers, template->atomNumbers, template->atomCount);
    memset(molecule->atomNumbers + template->atomCount, 1, template->hydrogenAtomCount);

    molecule->atomHydrogens = (uint8_t *) alloc_memory(&memory, molecule->atomCount);
    memcpy(molecule->atomHydrogens, template->atomHydrogens, template->atomCount);
    memset(molecule->atomHydrogens + template->atomCount, 0, template->hydrogenAtomCount);

    molecule->bondTypes = (uint8_t *) alloc_memory(&memory, molecule->bondCount);
    memcpy(molecule->bondTypes, template->bondTypes, template->bondCount);
    memset(molecule->bondTypes + template->bondCount, BOND_SINGLE, template->hydrogenAtomCount);

    if(template->atomCharges)
    {
        molecule->atomCharges = (int8_t *) alloc_memory(&memory, molecule->atomCount);
        memcpy(molecule->atomCharges, template->atomCharges, template->atomCount);
        memset(molecule->atomCharges + template->atomCount, 0, template->hydrogenAtomCount);
    }
    else
    {
        molecule->atomCharges = NULL;
    }

    if(template->atomMasses)
    {
        molecule->atomMasses = (int8_t *) alloc_memory(&memory, molecule->atomCount);
        memcpy(molecule->atomMasses, template->atomMasses, template->atomCount);
        memset(molecule->atomMasses + template->atomCount, 0, template->hydrogenAtomCount);
    }
    else
    {
        molecule->atomMasses = NULL;
    }

    if(template->atomRadicalTypes)
    {
        molecule->atomRadicalTypes = (int8_t *) alloc_memory(&memory, molecule->atomCount);
        memcpy(molecule->atomRadicalTypes, template->atomRadicalTypes, template->atomCount);
        memset(molecule->atomRadicalTypes + template->atomCount, 0, template->hydrogenAtomCount);
    }
    else
    {
        molecule->atomRadicalTypes = NULL;
    }

    if(template->atomStereo)
    {
        molecule->atomStereo = (uint8_t *) alloc_memory(&memory, molecule->atomCount);
        memcpy(molecule->atomStereo, template->atomStereo, template->atomCount);
        memset(molecule->atomStereo + template->atomCount, 0, template->hydrogenAtomCount);
    }
    else
    {
        molecule->atomStereo = NULL;
    }

    if(template->bondStereo)
    {
        molecule->bondStereo = (uint8_t *) alloc_memory(&memory, molecule->bondCount);
        memcpy(molecule->bondStereo, template->bondStereo, template->bondCount);
        memset(molecule->bondStereo + template->bondCount, 0, template->hydrogenAtomCount);
    }
    else
    {
        molecule->bondStereo = NULL;
    }

    molecule->bondLists = (AtomIdx *) alloc_memory(&memory, BOND_LIST_BASE_SIZE * molecule->atomCount * sizeof(AtomIdx));
    memcpy(molecule->bondLists, template->bondLists, BOND_LIST_BASE_SIZE * template->atomCount * sizeof(AtomIdx));

    molecule->bondListSizes = (MolSize *) alloc_memory_zero(&memory, molecule->atomCount * sizeof(MolSize));
    memcpy(molecule->bondListSizes, template->bondListSizes, template->atomCount * sizeof(MolSize));
    memset(molecule->bondListSizes + template->atomCount, 0, template->hydrogenAtomCount * sizeof(MolSize));

    molecule->contains = (AtomIdx (*)[2]) alloc_memory(&memory, molecule->bondCount * 2 * sizeof(AtomIdx));
    memcpy(molecule->contains, template->contains, template->bondCount * 2 * sizeof(AtomIdx));

    molecule->bondMatrix = (BondIdx *) alloc_memory_one(&memory, molecule->atomCount * molecule->atomCount * sizeof(BondIdx));

    for(int i = 0; i < template->atomCount; i++)
        for(int j = 0; j < template->atomCount; j++)
            molecule->bondMatrix[i * molecule->atomCount + j] = template->bondMatrix[i * template->atomCount + j];


    int hydrogen = template->atomCount;
    BondIdx boundIdx = template->bondCount;

    for(AtomIdx a = 0; a < template->atomCount; a++)
    {
        for(int h = 0; h < template->atomHydrogens[a]; h++)
        {
            molecule->bondLists[a * BOND_LIST_BASE_SIZE + molecule->bondListSizes[a]++] = hydrogen;
            molecule->bondLists[hydrogen * BOND_LIST_BASE_SIZE + molecule->bondListSizes[hydrogen]++] = a;

            if(unlikely(molecule->bondListSizes[hydrogen] == BOND_LIST_BASE_SIZE || molecule->bondListSizes[hydrogen] == BOND_LIST_BASE_SIZE))
                return NULL;

            molecule->bondMatrix[a * molecule->atomCount + hydrogen] = boundIdx;
            molecule->bondMatrix[hydrogen * molecule->atomCount + a] = boundIdx;

            molecule->contains[boundIdx][0] = a;
            molecule->contains[boundIdx][1] = hydrogen;

            boundIdx++;
            hydrogen++;
        }
    }


    molecule->labelCount = template->labelCount;
    molecule->labels = NULL;

    if(template->labels != NULL)
    {
        molecule->labels = (AtomLabel *) alloc_memory(&memory, template->labelCount * sizeof(AtomLabel));

        for(int i = 0; i < molecule->labelCount; i++)
        {
            molecule->labels[i].atom = template->labels[i].atom;
            molecule->labels[i].length = template->labels[i].length;
            molecule->labels[i].label = (uint8_t *) alloc_memory(&memory, template->labels[i].length);
            memcpy(molecule->labels[i].label, template->labels[i].label, template->labels[i].length);
        }
    }


    molecule->sgroupCount = template->sgroupCount;
    molecule->sgroups = NULL;

    if(template->sgroups != NULL)
    {
        molecule->sgroups = (SGroup *) alloc_memory(&memory, template->sgroupCount * sizeof(SGroup));

        for(int i = 0; i < molecule->sgroupCount; i++)
        {
            molecule->sgroups[i].type = template->sgroups[i].type;
            molecule->sgroups[i].subtype = template->sgroups[i].subtype;
            molecule->sgroups[i].connectivity = template->sgroups[i].connectivity;

            molecule->sgroups[i].atomCount = template->sgroups[i].atomCount;
            molecule->sgroups[i].atoms = (AtomIdx *) alloc_memory(&memory, template->sgroups[i].atomCount * sizeof(AtomIdx));
            memcpy(molecule->sgroups[i].atoms, template->sgroups[i].atoms, template->sgroups[i].atomCount * sizeof(AtomIdx));

            molecule->sgroups[i].bondCount = template->sgroups[i].bondCount;
            molecule->sgroups[i].bonds = (AtomIdx (*)[2]) alloc_memory(&memory, template->sgroups[i].bondCount * 2 * sizeof(AtomIdx));
            memcpy(molecule->sgroups[i].bonds, template->sgroups[i].bonds, template->sgroups[i].bondCount * 2 * sizeof(AtomIdx));
        }
    }

    return molecule;
}


static inline bool molecule_is_extended_search_needed(const uint8_t *restrict data, bool withRGroups,
        bool withCharges, bool withIsotopes, bool withRadicals)
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
        if(withRGroups && ((int8_t) data[i]) < 0)
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

        if(value == 0 || (value & 0xFFF) >= heavyAtomCount)
            return true;

        hBonds[i]++;
    }

    data += hAtomCount * HBOND_BLOCK_SIZE;


    for(int i = 0; i < hAtomCount; i++)
        if(hBonds[i] != 1)
            return true;

    if(!withCharges && !withIsotopes && !withRadicals)
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

            case RECORD_RADICAL:
                if(withRadicals && idx >= heavyAtomCount)
                    return true;
                break;
        }
    }

    return false;
}


static inline bool molecule_has_pseudo_atom(const uint8_t *restrict data)
{
    int xAtomCount = *data << 8 | *(data + 1);
    data += 10;

    for(int i = 0; i < xAtomCount; i++)
        if(((int8_t) data[i]) < 0)
            return true;

    return false;
}


static inline bool molecule_has_charged_hydrogen(const uint8_t *restrict data)
{
    int xAtomCount = data[0] << 8 | data[1];
    int cAtomCount = data[2] << 8 | data[3];
    int hAtomCount = data[4] << 8 | data[5];
    int xBondCount = data[6] << 8 | data[7];
    int specialCount = data[8] << 8 | data[9];

    int heavyAtomCount = xAtomCount + cAtomCount;
    int base = 10 + xAtomCount + xBondCount * BOND_BLOCK_SIZE + hAtomCount * HBOND_BLOCK_SIZE;

    for(int i = 0; i < specialCount; i++)
    {
        int offset = base + i * SPECIAL_BLOCK_SIZE;
        int value = data[offset + 0] * 256 | data[offset + 1];
        int idx = value & 0xFFF;

        if(data[offset] >> 4 == RECORD_CHARGE && idx >= heavyAtomCount)
            return true;
    }

    return false;
}


static inline bool molecule_has_hydrogen_isotope(const uint8_t *restrict data)
{
    int xAtomCount = data[0] << 8 | data[1];
    int cAtomCount = data[2] << 8 | data[3];
    int hAtomCount = data[4] << 8 | data[5];
    int xBondCount = data[6] << 8 | data[7];
    int specialCount = data[8] << 8 | data[9];

    int heavyAtomCount = xAtomCount + cAtomCount;
    int base = 10 + xAtomCount + xBondCount * BOND_BLOCK_SIZE + hAtomCount * HBOND_BLOCK_SIZE;

    for(int i = 0; i < specialCount; i++)
    {
        int offset = base + i * SPECIAL_BLOCK_SIZE;
        int value = data[offset + 0] * 256 | data[offset + 1];
        int idx = value & 0xFFF;

        if(data[offset] >> 4 == RECORD_ISOTOPE && idx >= heavyAtomCount)
            return true;
    }

    return false;
}


static inline bool molecule_has_hydrogen_radical(const uint8_t *restrict data)
{
    int xAtomCount = data[0] << 8 | data[1];
    int cAtomCount = data[2] << 8 | data[3];
    int hAtomCount = data[4] << 8 | data[5];
    int xBondCount = data[6] << 8 | data[7];
    int specialCount = data[8] << 8 | data[9];

    int heavyAtomCount = xAtomCount + cAtomCount;
    int base = 10 + xAtomCount + xBondCount * BOND_BLOCK_SIZE + hAtomCount * HBOND_BLOCK_SIZE;

    for(int i = 0; i < specialCount; i++)
    {
        int offset = base + i * SPECIAL_BLOCK_SIZE;
        int value = data[offset + 0] * 256 | data[offset + 1];
        int idx = value & 0xFFF;

        if(data[offset] >> 4 == RECORD_RADICAL && idx >= heavyAtomCount)
            return true;
    }

    return false;
}


static inline bool molecule_has_sgroup(const uint8_t *restrict data)
{
    int xAtomCount = data[0] << 8 | data[1];
    int hAtomCount = data[4] << 8 | data[5];
    int xBondCount = data[6] << 8 | data[7];
    int specialCount = data[8] << 8 | data[9];

    int base = 10 + xAtomCount + xBondCount * BOND_BLOCK_SIZE + hAtomCount * HBOND_BLOCK_SIZE;

    for(int i = 0; i < specialCount; i++)
    {
        int offset = base + i * SPECIAL_BLOCK_SIZE;

        if(data[offset] >> 4 == RECORD_SGROUP_COUNT)
            return true;
    }

    return false;
}


static inline bool molecule_has_multivalent_hydrogen(const uint8_t *restrict data)
{
    int xAtomCount = data[0] << 8 | data[1];
    int cAtomCount = data[2] << 8 | data[3];
    int xBondCount = data[6] << 8 | data[7];

    if(xBondCount == 0)
        return false;

    int offset = 10 + xAtomCount + (xBondCount - 1) * BOND_BLOCK_SIZE;

    int b0 = data[offset + 0];
    int b1 = data[offset + 1];
    int b2 = data[offset + 2];

    int x = b0 | (b1 << 4 & 0xF00);
    int y = b2 | (b1 << 8 & 0xF00);

    return x >= xAtomCount + cAtomCount || y >= xAtomCount + cAtomCount;
}


static inline bool molecule_is_pseudo_atom(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->atomNumbers[atom] < 0;
}


static inline int8_t molecule_get_atom_number(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->atomNumbers[atom];
}


static inline AtomLabel *molecule_get_atom_label(const Molecule *restrict molecule, AtomIdx atom)
{
    for(int i = 0; i < molecule->labelCount; i++)
        if(molecule->labels[i].atom == atom)
            return molecule->labels + i;

    return NULL;
}


static inline bool molecule_is_metal(const Molecule *restrict molecule, AtomIdx atom)
{
    int8_t number = molecule_get_atom_number(molecule, atom);

    return  (number > 2 && number < 5) || (number > 10 && number < 14) ||
            (number > 18 && number < 32) || (number > 36 && number < 51) ||
            (number > 54 && number < 85) || number > 86;
}


static inline bool molecule_is_halogen(const Molecule *restrict molecule, AtomIdx atom)
{
    int8_t number = molecule_get_atom_number(molecule, atom);

    return number == 9 || number == 17 || number == 35 || number == 53 || number == 85;
}


static inline uint8_t molecule_get_hydrogen_count(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->atomHydrogens[atom];
}


static inline int8_t molecule_get_formal_charge(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->atomCharges[atom];
}


static inline int8_t molecule_get_atom_mass(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->atomMasses[atom];
}


static inline int8_t molecule_get_atom_radical_type(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->atomRadicalTypes[atom];
}


static inline AtomIdx *molecule_get_bonded_atom_list(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->bondLists + atom * BOND_LIST_BASE_SIZE;
}


static inline MolSize molecule_get_bonded_atom_list_size(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->bondListSizes[atom];
}


static inline bool molecule_has_restH_flags(const Molecule *restrict molecule)
{
    return molecule->restH != NULL;
}


static inline bool molecule_get_atom_restH_flag(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->restH[atom];
}


static inline AtomIdx *molecule_bond_atoms(const Molecule *restrict molecule, BondIdx bond)
{
    return molecule->contains[bond];
}


static inline bool molecule_bond_contains(const Molecule *restrict molecule, BondIdx bond, AtomIdx atom)
{
    return molecule->contains[bond][0] == atom || molecule->contains[bond][1] == atom;
}


static inline AtomIdx molecule_get_other_bond_atom(const Molecule *restrict molecule, BondIdx bond, AtomIdx atom)
{
    if(molecule->contains[bond][0] == atom)
        return molecule->contains[bond][1];
    else if(molecule->contains[bond][1] == atom)
        return molecule->contains[bond][0];
    return -1;
}


static inline AtomIdx molecule_get_opposite_atom(const Molecule *restrict molecule, AtomIdx centre, AtomIdx atom)
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


static inline BondIdx molecule_get_bond(const Molecule *restrict molecule, AtomIdx i, AtomIdx j)
{
    return molecule->bondMatrix[i * molecule->atomCount + j];
}


static inline uint8_t molecule_get_bond_type(const Molecule *restrict molecule, BondIdx bond)
{
    return molecule->bondTypes[bond];
}


static inline uint8_t molecule_get_atom_stereo(const Molecule *restrict molecule, AtomIdx atom)
{
    return molecule->atomStereo[atom];
}


static inline uint8_t molecule_get_bond_stereo(const Molecule *restrict molecule, BondIdx bond)
{
    return molecule->bondStereo[bond];
}


static inline AtomIdx molecule_get_last_stereo_bond_ligand(const Molecule *restrict molecule, AtomIdx atom, AtomIdx other, AtomIdx ligand)
{
    MolSize listSize = molecule_get_bonded_atom_list_size(molecule, atom);

    if(listSize == 2)
    {
        return MAX_ATOM_IDX;
    }
    else if(listSize == 3)
    {
        for(int i = 0; i < 3; i++)
        {
            AtomIdx a = molecule_get_bonded_atom_list(molecule, atom)[i];

            if(a != other && a != ligand)
                return a;
        }
    }

    // unexpected branch
    return MAX_ATOM_IDX;
}


static inline AtomIdx molecule_get_last_chiral_ligand(const Molecule *restrict molecule, AtomIdx centre, AtomIdx *ligands)
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

    // unexpected branch
    return MAX_ATOM_IDX;
}


static inline uint8_t normalize_atom_stereo(AtomIdx indexes[4], uint8_t stereo)
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
        return ~stereo & 0x03;

    return stereo;
}


static inline uint8_t normalize_bond_stereo(AtomIdx indexes[4], uint8_t conformation)
{
    bool reverse = false;

    if(indexes[0] > indexes[1])
        reverse = !reverse;

    if(indexes[2] > indexes[3])
        reverse = !reverse;

    if(reverse)
        return ~conformation & 0x03;

    return conformation;
}

#endif /* MOLECULE_H_ */
