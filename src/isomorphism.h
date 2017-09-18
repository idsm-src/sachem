#ifndef ISOMORPHISM_H__
#define ISOMORPHISM_H__

#include <postgres.h>
#include <stdbool.h>
#include "molecule.h"

#define UNDEFINED_CORE              -1
#define MASKED_TARGET               -1
#define is_core_defined(value)      ((value) >= 0)
#define is_target_masked(value)     ((value) < 0)


typedef struct
{
    int target_range_len;
    int target_selector;
    int target_idx;
} VF2Undo;


typedef struct
{
    bool strictStereo;
    bool exact;

    const Molecule *restrict query;
    const Molecule *restrict target;

    int queryAtomCount;
    int targetAtomCount;

    int *restrict core_query;
    int *restrict core_target;
    int core_len;

    bool *restrict query_range_flags;
    bool *restrict query_range_flag;
    int *restrict query_order;
    int *restrict query_range_lens;
    int query_range_len;
    int query_idx;

    bool *restrict target_range_flag;
    int *restrict target_range_stack;
    int target_range_len;
    int target_selector;
    int target_idx;

    VF2Undo *undos;

    bool *restrict empty_flags;
} VF2State;


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


inline void vf2state_init(VF2State *const restrict vf2state, const Molecule *const restrict query, bool strictStereo, bool exact)
{
    int queryAtomCount = query->atomCount;

    vf2state->strictStereo = strictStereo;
    vf2state->exact = exact;
    vf2state->query = query;
    vf2state->queryAtomCount = queryAtomCount;
    vf2state->core_len = 0;
    vf2state->query_order = palloc(queryAtomCount * sizeof(int));
    vf2state->query_range_lens = palloc(queryAtomCount * sizeof(int));
    vf2state->query_range_len = 0;
    vf2state->undos = palloc(queryAtomCount * sizeof(VF2Undo));


    vf2state->core_query = palloc(queryAtomCount * sizeof(int));

    for(int i = 0; i < queryAtomCount; i++)
        vf2state->core_query[i] = UNDEFINED_CORE;


    vf2state->query_range_flags = palloc((queryAtomCount + 1) * queryAtomCount * sizeof(bool));
    vf2state->query_range_flag = vf2state->query_range_flags;

    for(int i = 0; i < queryAtomCount; i++)
        vf2state->query_range_flag[i] = 0;

    for(int core = 0; core < queryAtomCount; core++)
    {
        memcpy(vf2state->query_range_flag + queryAtomCount, vf2state->query_range_flag, sizeof(bool) * queryAtomCount);
        vf2state->query_range_flag += queryAtomCount;

        int idx = 0;

        if(vf2state->query_range_len > core)
        {
            while(idx < queryAtomCount && (is_core_defined(vf2state->core_query[idx]) || !vf2state->query_range_flag[idx]))
                idx++;
        }
        else
        {
            while(idx < queryAtomCount && is_core_defined(vf2state->core_query[idx]))
                idx++;
        }

        vf2state->query_order[core] = idx;
        vf2state->core_query[idx] = 0;

        if(!vf2state->query_range_flag[idx])
        {
            vf2state->query_range_flag[idx] = true;
            vf2state->query_range_len++;
        }


        int *restrict queryBondList = molecule_get_bond_list(query, idx);
        int queryBondListSize = molecule_get_bond_list_size(query, idx);

        for(int i = 0; i < queryBondListSize; i++)
        {
            int other = queryBondList[i];

            if(!vf2state->query_range_flag[other])
            {
                vf2state->query_range_flag[other] = true;
                vf2state->query_range_len++;
            }
        }

        vf2state->query_range_lens[core] = vf2state->query_range_len;
    }


    if(unlikely(exact))
    {
        vf2state->target_range_flag = palloc(queryAtomCount * sizeof(bool));
        vf2state->target_range_stack = palloc(queryAtomCount * sizeof(int));
        vf2state->core_target = palloc(queryAtomCount * sizeof(int));
    }
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
    if(vf2state->query_range_len > vf2state->core_len)
    {
        vf2state->target_selector++;

        while(vf2state->target_selector < vf2state->target_range_len)
        {
            if(!is_target_masked(vf2state->target_range_stack[vf2state->target_selector]))
            {
                vf2state->target_idx = vf2state->target_range_stack[vf2state->target_selector];
                return true;
            }

            vf2state->target_selector++;
        }
    }
    else
    {
        vf2state->target_idx++;

        while(vf2state->target_idx < vf2state->targetAtomCount)
        {
            if(!is_core_defined(vf2state->core_target[vf2state->target_idx]))
            {
                vf2state->target_selector = vf2state->target_range_len;

                if(vf2state->target_range_flag[vf2state->target_idx])
                {
                    for(int i = 0; i < vf2state->target_range_len; i++)
                    {
                        if(vf2state->target_range_stack[i] == vf2state->target_idx)
                        {
                            vf2state->target_selector = i;
                            break;
                        }
                    }
                }

                return true;
            }

            vf2state->target_idx++;
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
        if(unlikely(molecule_is_pseudo_atom(vf2state->target, targetAtom)))
            return true;
        else if(queryAtomNumber != Q_ATOM_NUMBER) // "Q" is a heteroatom (i.e. any atom except C or H)
            return true;
        else if(targetAtomNumber == C_ATOM_NUMBER || targetAtomNumber == H_ATOM_NUMBER)
            return false;
        else
            return true;
    }

    if(unlikely(molecule_is_pseudo_atom(vf2state->target, targetAtom)))
    {
        if(unlikely(molecule_is_pseudo_atom(vf2state->query, queryAtom)))
            return true;
        else if(targetAtomNumber != Q_ATOM_NUMBER) // "Q" is a heteroatom (i.e. any atom except C or H)
            return true;
        else if(queryAtomNumber == C_ATOM_NUMBER || queryAtomNumber == H_ATOM_NUMBER)
            return false;
        else
            return true;
    }

    return false;
}


inline bool vf2state_bond_matches(const VF2State *const restrict vf2state, int qIdx1, int qIdx2, int tIdx1, int tIdx2)
{
    int queryBond = molecule_get_bond(vf2state->query, qIdx1, qIdx2);
    int targetbond = molecule_get_bond(vf2state->target, tIdx1, tIdx2);

    if(likely(targetbond < 0))
        return false;

    if(likely(molecule_get_bond_data(vf2state->query, queryBond) != molecule_get_bond_data(vf2state->target, targetbond)))
        return false;

    return true;
}


inline bool vf2state_is_feasible_pair(const VF2State *const restrict vf2state)
{
    if(likely(!vf2state_atom_matches(vf2state, vf2state->query_idx, vf2state->target_idx)))
        return false;

    if(likely(!vf2state->exact))
    {
        if(unlikely(molecule_get_formal_charge(vf2state->query, vf2state->query_idx) != 0 &&
                molecule_get_formal_charge(vf2state->query, vf2state->query_idx) !=
                molecule_get_formal_charge(vf2state->target, vf2state->target_idx)))
            return false;

        if(unlikely(molecule_get_hydrogen_count(vf2state->query, vf2state->query_idx) >
                molecule_get_hydrogen_count(vf2state->target, vf2state->target_idx) &&
                !molecule_is_pseudo_atom(vf2state->query, vf2state->query_idx) &&
                !molecule_is_pseudo_atom(vf2state->target, vf2state->target_idx)))
            return false;
    }
    else
    {
        if(unlikely(molecule_get_formal_charge(vf2state->query, vf2state->query_idx) !=
                molecule_get_formal_charge(vf2state->target, vf2state->target_idx)))
            return false;

        if(unlikely(molecule_get_hydrogen_count(vf2state->query, vf2state->query_idx) !=
                molecule_get_hydrogen_count(vf2state->target, vf2state->target_idx) &&
                !molecule_is_pseudo_atom(vf2state->query, vf2state->query_idx) &&
                !molecule_is_pseudo_atom(vf2state->target, vf2state->target_idx)))
            return false;
    }


    int terminQuery = 0;
    int newQuery = 0;
    int terminTarget = 0;
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
            if(vf2state->query_range_flag[other1])
                terminQuery++;
            else
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
            if(unlikely(vf2state->exact))
            {
                int other1 = vf2state->core_target[other2];

                if(!vf2state_bond_matches(vf2state, vf2state->query_idx, other1, vf2state->target_idx, other2))
                    return false;
            }
        }
        else
        {
            if(vf2state->target_range_flag[other2])
                terminTarget++;
            else
                newTarget++;
        }
    }

    if(unlikely(vf2state->exact))
        return terminQuery == terminTarget && newQuery == newTarget;
    else
        return terminQuery <= terminTarget && terminQuery + newQuery <= terminTarget + newTarget;
}


inline void vf2state_undo_add_pair(VF2State *const restrict vf2state)
{
    VF2Undo *restrict undo = &vf2state->undos[--vf2state->core_len];

    vf2state->query_idx = vf2state->query_order[vf2state->core_len];

    vf2state->core_query[vf2state->query_idx] = UNDEFINED_CORE;
    vf2state->core_target[undo->target_idx] = UNDEFINED_CORE;

    vf2state->target_range_stack[undo->target_selector] = undo->target_idx;

    for(int i = undo->target_range_len; i < vf2state->target_range_len; i++)
    {
        if(is_core_defined(vf2state->target_range_stack[i]))
            vf2state->target_range_flag[vf2state->target_range_stack[i]] = false;
        else
            vf2state->target_range_flag[undo->target_idx] = false;
    }

    vf2state->target_range_len = undo->target_range_len;
    vf2state->target_selector = undo->target_selector;
    vf2state->target_idx = undo->target_idx;

    vf2state->query_range_len = vf2state->core_len > 0 ? vf2state->query_range_lens[vf2state->core_len - 1] : 0;
    vf2state->query_range_flag = vf2state->query_range_flags + vf2state->queryAtomCount * vf2state->core_len;
}


inline void vf2state_add_pair(VF2State *const restrict vf2state)
{
    VF2Undo *restrict undo = &vf2state->undos[vf2state->core_len];

    vf2state->query_range_len = vf2state->query_range_lens[vf2state->core_len];
    vf2state->query_range_flag = vf2state->query_range_flags + vf2state->queryAtomCount * (vf2state->core_len + 1);

    vf2state->core_len++;
    vf2state->core_query[vf2state->query_idx] = vf2state->target_idx;
    vf2state->core_target[vf2state->target_idx] = vf2state->query_idx;

    undo->target_range_len = vf2state->target_range_len;
    undo->target_selector = vf2state->target_selector;
    undo->target_idx = vf2state->target_idx;

    if(!vf2state->target_range_flag[vf2state->target_idx])
    {
        vf2state->target_range_flag[vf2state->target_idx] = true;
        vf2state->target_selector = vf2state->target_range_len++;
    }

    vf2state->target_range_stack[vf2state->target_selector] = MASKED_TARGET;


    int *restrict targetBondList = molecule_get_bond_list(vf2state->target, vf2state->target_idx);
    int targetBondListSize = molecule_get_bond_list_size(vf2state->target, vf2state->target_idx);

    for(int i = 0; i < targetBondListSize; i++)
    {
        int other = targetBondList[i];

        if(!vf2state->target_range_flag[other])
        {
            vf2state->target_range_flag[other] = true;
            vf2state->target_range_stack[vf2state->target_range_len++] = other;
        }
    }
}


bool vf2state_is_stereo_valid(const VF2State *const restrict vf2state)
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

        if(queryStereo != NONE_STEREO && queryStereo != UNDEF_STEREO)
        {
            if(targetStereo == NONE_STEREO)
                continue;

            if(targetStereo == UNDEF_STEREO)
                return false;


            int queryAtoms[4];
            int listSize = molecule_get_bond_list_size(query, queryAtomIdx);

            for(int i = 0; i < listSize; i++)
                queryAtoms[i] = molecule_get_bond_list(query, queryAtomIdx)[i];

            if(listSize == 3)
                queryAtoms[3] = queryAtomIdx;

            sort_stereo_atoms(queryAtoms);


            int targetAtoms[4];

            for(int i = 0; i < 4; i++)
                targetAtoms[i] = vf2state->core_query[queryAtoms[i]];

            if(normalize_atom_stereo(targetAtoms, targetStereo) != queryStereo)
                return false;
        }
        else if(vf2state->exact)
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

        if(queryStereo != NONE_STEREO && queryStereo != UNDEF_STEREO)
        {
            if(targetStereo == NONE_STEREO)
                continue;

            if(targetStereo == UNDEF_STEREO)
                return false;


            int queryAtoms[4];

            int *bondList0 = molecule_get_bond_list(query, queryBondAtoms[0]);
            int listSize0 = molecule_get_bond_list_size(query, queryBondAtoms[0]);

            int idx = 0;

            for(int i = 0; i < listSize0; i++)
                if(bondList0[i] != queryBondAtoms[1])
                    queryAtoms[idx++] = bondList0[i];

            if(listSize0 == 2)
                queryAtoms[idx++] = queryBondAtoms[0];

            int *bondList1 = molecule_get_bond_list(query, queryBondAtoms[1]);
            int listSize1 = molecule_get_bond_list_size(query, queryBondAtoms[1]);

            for(int i = 0; i < listSize1; i++)
                if(bondList1[i] != queryBondAtoms[0])
                    queryAtoms[idx++] = bondList1[i];

            if(listSize1 == 2)
                queryAtoms[idx] = queryBondAtoms[1];

            sort_bond_atoms(queryAtoms);


            int targetAtoms[4];

            for(int i = 0; i < 4; i++)
                targetAtoms[i] = vf2state->core_query[queryAtoms[i]];

            if(normalize_bond_stereo(targetAtoms, targetStereo) != queryStereo)
                return false;
        }
        else if(vf2state->exact)
        {
            if(queryStereo != targetStereo)
                return false;
        }
    }


    return true;
}


bool vf2state_is_match_valid(const VF2State *const restrict vf2state)
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


    if(unlikely(vf2state->strictStereo) && !vf2state_is_stereo_valid(vf2state))
            return false;


    return true;
}


bool vf2state_match_core(VF2State *const restrict vf2state)
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

        if(vf2state->query_range_len > vf2state->target_range_len || !vf2state_next_query(vf2state))
            goto recursion_return;


        while(vf2state_next_target(vf2state))
        {
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


inline bool vf2state_match(VF2State *const restrict vf2state, const Molecule *const restrict target)
{
    if(likely(!vf2state->exact))
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
    vf2state->query_range_flag = vf2state->query_range_flags;
    vf2state->query_range_len = 0;
    vf2state->target_range_len = 0;

    if(likely(!vf2state->exact))
    {
        vf2state->target_range_flag = palloc(targetAtomCount * sizeof(bool));
        vf2state->target_range_stack = palloc(targetAtomCount * sizeof(int));
        vf2state->core_target = palloc(targetAtomCount * sizeof(int));
    }

    for(int i = 0; i < targetAtomCount; i++)
        vf2state->target_range_flag[i] = false;

    for(int i = 0; i < targetAtomCount; i++)
        vf2state->core_target[i] = UNDEFINED_CORE;

    for(int i = 0; i < vf2state->queryAtomCount; i++)
        vf2state->core_query[i] = UNDEFINED_CORE;


    return vf2state_match_core(vf2state);
}

#endif /* ISOMORPHISM_H__ */
