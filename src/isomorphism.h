/*
 * Copyright (C) 2001
 *   Dipartimento di Informatica e Sistemistica,
 *   Universita degli studi di Napoli ``Federico II'
 *   <http://amalfi.dis.unina.it>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 *  1. The above copyright notice and this permission notice shall be included in all copies or substantial
 *      portions of the Software, together with the associated disclaimers.
 *  2. Any modification to the standard distribution of the Software shall be mentioned in a prominent notice
 *      in the documentation provided with the modified distribution, stating clearly how, when and by
 *      whom the Software has been modified.
 *  3. Either the modified distribution shall contain the entire sourcecode of the standard distribution of the
 *      Software, or the documentation shall provide instructions on where the source code of the standard
 *      distribution of the Software can be obtained.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef ISOMORPHISM_H__
#define ISOMORPHISM_H__

#include <postgres.h>
#include <utils/timeout.h>
#include <stdbool.h>
#include "molecule.h"

#define USE_VF2_TIMEOUT              1

#define UNDEFINED_CORE              -1
#define MASKED_TARGET               -1
#define is_core_defined(value)      ((value) >= 0)
#define is_target_masked(value)     ((value) < 0)


typedef enum
{
    GRAPH_SUBSTRUCTURE = 0,
    GRAPH_EXACT = 1
} GraphMode;


typedef enum
{
    CHARGE_IGNORE = 0,
    CHARGE_DEFAULT_AS_UNCHARGED = 1,
    CHARGE_DEFAULT_AS_ANY = 2
} ChargeMode;


typedef enum
{
    ISOTOPE_IGNORE = 0,
    ISOTOPE_DEFAULT_AS_STANDARD = 1,
    ISOTOPE_DEFAULT_AS_ANY = 2
} IsotopeMode;


typedef enum
{
    STEREO_IGNORE = 0,
    STEREO_STRICT = 1
} StereoMode;


typedef struct
{
    int target_selector;
    int target_idx;
} VF2Undo;


typedef struct
{
    GraphMode graphMode;
    ChargeMode chargeMode;
    IsotopeMode isotopeMode;
    StereoMode stereoMode;

    int targetId;

    const Molecule *restrict query;
    const Molecule *restrict target;

    int queryAtomCount;
    int targetAtomCount;

    int *restrict core_query;
    int *restrict core_target;
    int core_len;

    int *restrict query_order;
    int *restrict query_parents;
    int query_idx;

    int target_selector;
    int target_idx;

    VF2Undo *undos;
} VF2State;


#if USE_VF2_TIMEOUT
extern TimeoutId vf2TimeoutId;
extern volatile bool vf2Timeouted;

void vf2_timeout_handler(void);
#endif


inline void swap_int(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}


inline void sort_stereo_atoms(int array[4])
{
    for(int n = 4; n > 0; n--)
        for(int i = 1; i < n; i++)
            if(array[i-1] > array[i])
                swap_int(array + i, array + i - 1);
}


inline void sort_bond_atoms(int array[4])
{
    if(array[0] > array[1])
        swap_int(array + 0, array + 1);

    if(array[2] > array[3])
        swap_int(array + 2, array + 3);
}


inline void vf2state_init(VF2State *const restrict vf2state, const Molecule *const restrict query,
        GraphMode graphMode, ChargeMode chargeMode, IsotopeMode isotopeMode, StereoMode stereoMode)
{
#if USE_VF2_TIMEOUT
    if(unlikely(vf2TimeoutId == MAX_TIMEOUTS))
        vf2TimeoutId = RegisterTimeout(USER_TIMEOUT, vf2_timeout_handler);
#endif

    int queryAtomCount = query->atomCount;

    vf2state->graphMode = graphMode;
    vf2state->chargeMode = chargeMode;
    vf2state->isotopeMode = isotopeMode;
    vf2state->stereoMode = stereoMode;
    vf2state->query = query;
    vf2state->queryAtomCount = queryAtomCount;
    vf2state->core_len = 0;
    vf2state->core_query = (int *) palloc(queryAtomCount * sizeof(int));
    vf2state->query_order = (int *) palloc(queryAtomCount * sizeof(int));
    vf2state->query_parents = (int *) palloc(queryAtomCount * sizeof(int));
    vf2state->undos = (VF2Undo *) palloc(queryAtomCount * sizeof(VF2Undo));


    for(int i = 0; i < queryAtomCount; i++)
        vf2state->query_parents[i] = -1;

    uint8_t query_flags[queryAtomCount];

    for(int i = 0; i < queryAtomCount; i++)
        query_flags[i] = 0;

    for(int idx = 0; idx < queryAtomCount; idx++)
    {
        int selected = -1;
        int fallback = -1;

        for(int i = 0; i < queryAtomCount; i++)
        {
            if(selected == -1 && query_flags[i] == 1)
            {
                selected = i;
                break;
            }

            if(fallback == -1 && query_flags[i] == 0)
                fallback = i;
        }

        if(selected == -1)
            selected = fallback;


        int *restrict queryBondList = molecule_get_bond_list(query, selected);
        int queryBondListSize = molecule_get_bond_list_size(query, selected);

        query_flags[selected] = 2;

        for(int i = 0; i < queryBondListSize; i++)
        {
            int idx = queryBondList[i];

            if(query_flags[idx] == 0)
            {
                query_flags[idx] = 1;
                vf2state->query_parents[idx] = selected;
            }
        }

        vf2state->query_order[idx] = selected;
    }


    if(unlikely(graphMode == GRAPH_EXACT))
        vf2state->core_target = (int *) palloc(queryAtomCount * sizeof(int));
}


inline bool vf2state_next_query(VF2State *const restrict vf2state)
{
    if(unlikely(vf2state->core_len >= vf2state->queryAtomCount))
        return false;

    vf2state->query_idx = vf2state->query_order[vf2state->core_len];
    vf2state->target_idx = -1;
    vf2state->target_selector = -1;

    return true;
}


inline bool vf2state_next_target(VF2State *const restrict vf2state)
{
    int query_parent = vf2state->query_parents[vf2state->query_idx];

    if(likely(query_parent >= 0))
    {
        int target_parent = vf2state->core_query[query_parent];
        int *restrict targetBondList = molecule_get_bond_list(vf2state->target, target_parent);
        int targetBondListSize = molecule_get_bond_list_size(vf2state->target, target_parent);

        for(vf2state->target_selector++; vf2state->target_selector < targetBondListSize; vf2state->target_selector++)
        {
            int target_id = targetBondList[vf2state->target_selector];

            if(!is_core_defined(vf2state->core_target[target_id]))
            {
                vf2state->target_idx = target_id;
                return true;
            }
        }
    }
    else
    {
        for(vf2state->target_idx++; vf2state->target_idx < vf2state->targetAtomCount; vf2state->target_idx++)
        {
            if(!is_core_defined(vf2state->core_target[vf2state->target_idx]))
                return true;
        }
    }

    return false;
}


inline bool vf2state_atom_matches(const VF2State *const restrict vf2state, int queryAtom, int targetAtom)
{
    int queryAtomNumber = molecule_get_atom_number(vf2state->query, queryAtom);
    int targetAtomNumber = molecule_get_atom_number(vf2state->target, targetAtom);

    if(unlikely(queryAtomNumber == targetAtomNumber))
        return true;

    if(unlikely(molecule_is_pseudo_atom(vf2state->query, queryAtom)))
    {
        switch(queryAtomNumber)
        {
            case Q_ATOM_NUMBER:
                return targetAtomNumber != C_ATOM_NUMBER && targetAtomNumber != H_ATOM_NUMBER;

            case M_ATOM_NUMBER:
                return molecule_is_metal(vf2state->target, targetAtom) ||
                        molecule_is_pseudo_atom(vf2state->target, targetAtom) && targetAtom != X_ATOM_NUMBER;

            case X_ATOM_NUMBER:
                return molecule_is_halogen(vf2state->target, targetAtom) ||
                        molecule_is_pseudo_atom(vf2state->target, targetAtom) && targetAtom != M_ATOM_NUMBER;

            default:
                return true;
        }
    }

    if(unlikely(molecule_is_pseudo_atom(vf2state->target, targetAtom)))
    {
        switch(targetAtomNumber)
        {
            case Q_ATOM_NUMBER:
                return queryAtomNumber != C_ATOM_NUMBER && queryAtomNumber != H_ATOM_NUMBER;

            case M_ATOM_NUMBER:
                return molecule_is_metal(vf2state->query, queryAtom) ||
                        molecule_is_pseudo_atom(vf2state->query, queryAtom) && queryAtom != X_ATOM_NUMBER;

            case X_ATOM_NUMBER:
                return molecule_is_halogen(vf2state->query, queryAtom) ||
                        molecule_is_pseudo_atom(vf2state->query, queryAtom) && queryAtom != M_ATOM_NUMBER;

            default:
                return true;
        }
    }

    return false;
}


inline bool vf2state_bond_matches(const VF2State *const restrict vf2state, int qIdx1, int qIdx2, int tIdx1, int tIdx2)
{
    int queryBond = molecule_get_bond(vf2state->query, qIdx1, qIdx2);
    int targetbond = molecule_get_bond(vf2state->target, tIdx1, tIdx2);

    if(likely(targetbond < 0))
        return false;

    int queryBondType = molecule_get_bond_type(vf2state->query, queryBond);
    int targetbondType = molecule_get_bond_type(vf2state->target, targetbond);

    if(likely(queryBondType <= BOND_AROMATIC && targetbondType <= BOND_AROMATIC))
        return queryBondType == targetbondType;

    if(queryBondType > BOND_AROMATIC && targetbondType > BOND_AROMATIC)
        return true;

    if(queryBondType == BOND_ANY)
        return true;
    else if(queryBondType == BOND_SINGLE_OR_DOUBLE)
        return targetbondType == BOND_SINGLE || targetbondType == BOND_DOUBLE;
    else if(queryBondType == BOND_SINGLE_OR_AROMATIC)
        return targetbondType == BOND_SINGLE || targetbondType == BOND_AROMATIC;
    else if(queryBondType == BOND_DOUBLE_OR_AROMATIC)
        return targetbondType == BOND_DOUBLE || targetbondType == BOND_AROMATIC;

    if(targetbondType == BOND_ANY)
        return true;
    else if(targetbondType == BOND_SINGLE_OR_DOUBLE)
        return queryBondType == BOND_SINGLE || queryBondType == BOND_DOUBLE;
    else if(targetbondType == BOND_SINGLE_OR_AROMATIC)
        return queryBondType == BOND_SINGLE || queryBondType == BOND_AROMATIC;
    else if(targetbondType == BOND_DOUBLE_OR_AROMATIC)
        return queryBondType == BOND_DOUBLE || queryBondType == BOND_AROMATIC;

    return false;
}


inline bool vf2state_is_feasible_pair(const VF2State *const restrict vf2state)
{
    if(likely(!vf2state_atom_matches(vf2state, vf2state->query_idx, vf2state->target_idx)))
        return false;


    if(vf2state->chargeMode != CHARGE_IGNORE)
    {
        int queryCharge = molecule_get_formal_charge(vf2state->query, vf2state->query_idx);
        int targetCharge = molecule_get_formal_charge(vf2state->target, vf2state->target_idx);

        if(queryCharge != targetCharge && (queryCharge != 0 || vf2state->chargeMode == CHARGE_DEFAULT_AS_UNCHARGED))
            return false;
    }


    if(vf2state->isotopeMode != ISOTOPE_IGNORE)
    {
        int queryMass = molecule_get_atom_mass(vf2state->query, vf2state->query_idx);
        int targetMass = molecule_get_atom_mass(vf2state->target, vf2state->target_idx);

        if(queryMass != targetMass && (queryMass != 0 || vf2state->isotopeMode == ISOTOPE_DEFAULT_AS_STANDARD))
            return false;
    }


    if(likely(vf2state->graphMode != GRAPH_EXACT))
    {
        if(unlikely(molecule_get_hydrogen_count(vf2state->query, vf2state->query_idx) >
                molecule_get_hydrogen_count(vf2state->target, vf2state->target_idx) &&
                !molecule_is_pseudo_atom(vf2state->query, vf2state->query_idx) &&
                !molecule_is_pseudo_atom(vf2state->target, vf2state->target_idx)))
            return false;
    }
    else
    {
        if(unlikely(molecule_get_hydrogen_count(vf2state->query, vf2state->query_idx) !=
                molecule_get_hydrogen_count(vf2state->target, vf2state->target_idx) &&
                !molecule_is_pseudo_atom(vf2state->query, vf2state->query_idx) &&
                !molecule_is_pseudo_atom(vf2state->target, vf2state->target_idx)))
            return false;
    }


    int newQuery = 0;
    int newTarget = 0;

    int *restrict queryBondList = molecule_get_bond_list(vf2state->query, vf2state->query_idx);
    int queryBondListSize = molecule_get_bond_list_size(vf2state->query, vf2state->query_idx);

    for(int i = 0; i < queryBondListSize; i++)
    {
        int other1 = queryBondList[i];

        if(is_core_defined(vf2state->core_query[other1]))
        {
            int other2 = vf2state->core_query[other1];

            if(!vf2state_bond_matches(vf2state, vf2state->query_idx, other1, vf2state->target_idx, other2))
                return false;
        }
        else
        {
            newQuery++;
        }
    }


    int *restrict targetBondList = molecule_get_bond_list(vf2state->target, vf2state->target_idx);
    int targetBondListSize = molecule_get_bond_list_size(vf2state->target, vf2state->target_idx);

    for(int i = 0; i < targetBondListSize; i++)
    {
        int other2 = targetBondList[i];

        if(is_core_defined(vf2state->core_target[other2]))
        {
            if(unlikely(vf2state->graphMode == GRAPH_EXACT))
            {
                int other1 = vf2state->core_target[other2];

                if(!vf2state_bond_matches(vf2state, vf2state->query_idx, other1, vf2state->target_idx, other2))
                    return false;
            }
        }
        else
        {
            newTarget++;
        }
    }

    if(unlikely(vf2state->graphMode == GRAPH_EXACT))
        return newQuery == newTarget;
    else
        return newQuery <= newTarget;
}


inline void vf2state_undo_add_pair(VF2State *const restrict vf2state)
{
    VF2Undo *restrict undo = &vf2state->undos[--vf2state->core_len];

    vf2state->query_idx = vf2state->query_order[vf2state->core_len];

    vf2state->core_query[vf2state->query_idx] = UNDEFINED_CORE;
    vf2state->core_target[undo->target_idx] = UNDEFINED_CORE;

    vf2state->target_selector = undo->target_selector;
    vf2state->target_idx = undo->target_idx;
}


inline void vf2state_add_pair(VF2State *const restrict vf2state)
{
    VF2Undo *restrict undo = &vf2state->undos[vf2state->core_len];

    vf2state->core_len++;
    vf2state->core_query[vf2state->query_idx] = vf2state->target_idx;
    vf2state->core_target[vf2state->target_idx] = vf2state->query_idx;

    undo->target_selector = vf2state->target_selector;
    undo->target_idx = vf2state->target_idx;
}


inline bool vf2state_is_stereo_valid(const VF2State *const restrict vf2state)
{
    const Molecule *const restrict query = vf2state->query;
    const Molecule *const restrict target = vf2state->target;

    int queryAtomCount = query->atomCount;
    int targetAtomCount = target->atomCount;

    int queryBondCount = query->bondCount;
    int targetBondCount = target->bondCount;


    for(int queryAtomIdx = 0; queryAtomIdx < queryAtomCount; queryAtomIdx++)
    {
        int queryStereo = molecule_get_atom_stereo(query, queryAtomIdx);

        int targetAtomIdx = vf2state->core_query[queryAtomIdx];
        int targetStereo = molecule_get_atom_stereo(target, targetAtomIdx);

        if(queryStereo != TETRAHEDRAL_STEREO_NONE && queryStereo != TETRAHEDRAL_STEREO_UNDEFINED)
        {
            if(targetStereo == TETRAHEDRAL_STEREO_NONE)
                continue;

            if(targetStereo == TETRAHEDRAL_STEREO_UNDEFINED)
                return false;


            int queryAtoms[4];
            int listSize = molecule_get_bond_list_size(query, queryAtomIdx);

            for(int i = 0; i < listSize; i++)
                queryAtoms[i] = molecule_get_bond_list(query, queryAtomIdx)[i];

            if(listSize == 3)
                queryAtoms[3] = INT_MAX;

            sort_stereo_atoms(queryAtoms);


            int targetAtoms[4];

            for(int i = 0; i < listSize; i++)
                targetAtoms[i] = vf2state->core_query[queryAtoms[i]];

            if(listSize == 3)
                targetAtoms[3] = molecule_get_last_chiral_ligand(target, vf2state->core_query[queryAtomIdx], targetAtoms);

            if(normalize_atom_stereo(targetAtoms, targetStereo) != queryStereo)
                return false;
        }
        else if(vf2state->graphMode == GRAPH_EXACT)
        {
            if(queryStereo != targetStereo)
                return false;
        }
    }


    for(int queryBondIdx = 0; queryBondIdx < queryBondCount; queryBondIdx++)
    {
        int queryStereo = molecule_get_bond_stereo(query, queryBondIdx);

        int *queryBondAtoms = molecule_bond_atoms(query, queryBondIdx);
        int targetBondAtom0 = vf2state->core_query[queryBondAtoms[0]];
        int targetBondAtom1 = vf2state->core_query[queryBondAtoms[1]];
        int targetBondIdx = molecule_get_bond(target, targetBondAtom0, targetBondAtom1);
        int targetStereo = molecule_get_bond_stereo(target, targetBondIdx);

        if(queryStereo != BOND_STEREO_NONE && queryStereo != BOND_STEREO_UNDEFINED)
        {
            if(targetStereo == BOND_STEREO_NONE)
                continue;

            if(targetStereo == BOND_STEREO_UNDEFINED)
                return false;


            int queryAtoms[4];

            int *bondList0 = molecule_get_bond_list(query, queryBondAtoms[0]);
            int listSize0 = molecule_get_bond_list_size(query, queryBondAtoms[0]);

            int idx = 0;

            for(int i = 0; i < listSize0; i++)
                if(bondList0[i] != queryBondAtoms[1])
                    queryAtoms[idx++] = bondList0[i];

            if(listSize0 == 2)
                queryAtoms[idx++] = INT_MAX;

            int *bondList1 = molecule_get_bond_list(query, queryBondAtoms[1]);
            int listSize1 = molecule_get_bond_list_size(query, queryBondAtoms[1]);

            for(int i = 0; i < listSize1; i++)
                if(bondList1[i] != queryBondAtoms[0])
                    queryAtoms[idx++] = bondList1[i];

            if(listSize1 == 2)
                queryAtoms[idx] =  INT_MAX;

            sort_bond_atoms(queryAtoms);


            int targetAtoms[4];

            for(int i = 0; i < 4; i++)
                if(queryAtoms[i] != INT_MAX)
                    targetAtoms[i] = vf2state->core_query[queryAtoms[i]];

            if(queryAtoms[1] == INT_MAX)
                targetAtoms[1] = molecule_get_last_stereo_bond_ligand(target, targetBondAtom0, targetBondIdx, targetAtoms[0]);

            if(queryAtoms[3] == INT_MAX)
                targetAtoms[3] = molecule_get_last_stereo_bond_ligand(target, targetBondAtom1, targetBondIdx, targetAtoms[2]);

            if(normalize_bond_stereo(targetAtoms, targetStereo) != queryStereo)
                return false;
        }
        else if(vf2state->graphMode == GRAPH_EXACT)
        {
            if(queryStereo != targetStereo)
                return false;
        }
    }


    return true;
}


inline bool vf2state_is_match_valid(const VF2State *const restrict vf2state)
{
    const Molecule *const restrict query = vf2state->query;
    const Molecule *const restrict target = vf2state->target;

    int queryAtomCount = query->atomCount;
    int targetAtomCount = target->atomCount;

    int queryBondCount = query->bondCount;
    int targetBondCount = target->bondCount;

    /*
    The goal has been reached, so the query has been mapped to target,
    therfore query is a substructure of the target.

    However, if this was an R-group query, the result could still be
    rejected.If the RestH property is true for some atom with an R-group
    linked, then the R-group may only be substituted with a member
    of the Rgroup or with H..
    This can be verified:
        - find any atom in the query with RestH flagged
        - find the atom mapped to it in the target container
        - see if the target has more (non hydrogen) bonds than the query.
          if so,discard it.
    */

    if(molecule_has_restH_flags(query))
    {
        for(int i = 0; i < queryAtomCount; i++)
        {
            if(molecule_get_atom_restH_flag(query, i) == true)
            {
                int targetAtomIdx = vf2state->core_query[i];

                int queryConnectivityCount = 0;
                int targetConnectivityCount = 0;

                for(int b = 0; b < queryBondCount; b++)
                    if(molecule_bond_contains(query, b, i) &&
                            !(molecule_get_atom_number(query, molecule_get_bond_connected_atom(query, b, i)) == H_ATOM_NUMBER))
                        queryConnectivityCount++;

                for(int b = 0; b < targetBondCount; b++)
                    if(molecule_bond_contains(target, b, targetAtomIdx) &&
                            !(molecule_get_atom_number(target, molecule_get_bond_connected_atom(target, b, targetAtomIdx)) == H_ATOM_NUMBER))
                        targetConnectivityCount++;

                if(targetConnectivityCount > queryConnectivityCount)
                    return false;
            }
        }
    }


    if(unlikely(vf2state->stereoMode == STEREO_STRICT) && !vf2state_is_stereo_valid(vf2state))
            return false;


    return true;
}


inline bool vf2state_match_core(VF2State *const restrict vf2state)
{
    while(true)
    {
        recursion_entry:

        if(unlikely(vf2state->core_len == vf2state->query->atomCount))
        {
            if(vf2state_is_match_valid(vf2state))
                return true;

            goto recursion_return;
        }

        if(!vf2state_next_query(vf2state))
            goto recursion_return;


        while(vf2state_next_target(vf2state))
        {
            if(unlikely(QueryCancelPending))
                return false;

#if USE_VF2_TIMEOUT
            if(unlikely(vf2Timeouted))
            {
                elog(WARNING, "isomorphism: VF2 timeout expired for target %i", vf2state->targetId);
                return false;
            }
#endif

            if(vf2state_is_feasible_pair(vf2state))
            {
                vf2state_add_pair(vf2state);
                goto recursion_entry;

                recursion_return:

                if(vf2state->core_len == 0)
                    return false;

                vf2state_undo_add_pair(vf2state);
            }
        }

        goto recursion_return;
    }

    return false;
}


inline bool vf2state_match(VF2State *const restrict vf2state, const Molecule *const restrict target, int targetId, int timeout)
{
    vf2state->targetId = targetId;

    if(likely(vf2state->graphMode != GRAPH_EXACT))
    {
        if(vf2state->queryAtomCount > target->atomCount || vf2state->query->bondCount > target->bondCount)
            return false;
    }
    else
    {
        if(vf2state->queryAtomCount != target->atomCount || vf2state->query->bondCount != target->bondCount)
            return false;
    }

    int targetAtomCount = target->atomCount;

    vf2state->target = target;
    vf2state->targetAtomCount = targetAtomCount;
    vf2state->core_len = 0;

    if(likely(vf2state->graphMode != GRAPH_EXACT))
        vf2state->core_target = (int *) palloc(targetAtomCount * sizeof(int));

    for(int i = 0; i < targetAtomCount; i++)
        vf2state->core_target[i] = UNDEFINED_CORE;

    for(int i = 0; i < vf2state->queryAtomCount; i++)
        vf2state->core_query[i] = UNDEFINED_CORE;


#if USE_VF2_TIMEOUT
    vf2Timeouted = false;
    enable_timeout_after(vf2TimeoutId, timeout);
#endif
    bool result = vf2state_match_core(vf2state);
#if USE_VF2_TIMEOUT
    disable_timeout(vf2TimeoutId, false);
#endif

    return result;
}

#endif /* ISOMORPHISM_H__ */
