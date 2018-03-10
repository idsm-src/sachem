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
#ifndef SUBSTRUCTURE_MATCH_HPP__
#define SUBSTRUCTURE_MATCH_HPP__

#include <vector>
#include <algorithm>

extern "C"
{
#include "molecule.h"
}


#define UNDEFINED_CORE              -1
#define MASKED_TARGET               -1
#define is_core_defined(value)      ((value) >= 0)
#define is_target_masked(value)     ((value) < 0)


class SubstructureMatch
{
    private:
        class VF2Undo
        {
            public:
                int target_selector;
                int target_idx;
        };


        const Molecule *query;
        const Molecule *target;

        int queryAtomCount;
        int targetAtomCount;

        std::vector<int> core_query;
        std::vector<int> core_target;
        int core_len;

        std::vector<int> query_order;
        std::vector<int> query_parents;
        int query_idx;

        int target_selector;
        int target_idx;

        std::vector<VF2Undo> undos;
        std::vector<std::vector<int>> matches;
        size_t limit;


    public:
        SubstructureMatch(const Molecule *const query) :
            query(query),
            queryAtomCount(query->atomCount),
            core_query(query->atomCount),
            core_len(0),
            query_order(query->atomCount),
            query_parents(query->atomCount),
            undos(query->atomCount)
        {
            for(int i = 0; i < queryAtomCount; i++)
                query_parents[i] = -1;

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


                AtomIdx *queryBondList = molecule_get_bond_list(query, selected);
                MolSize queryBondListSize = molecule_get_bond_list_size(query, selected);

                query_flags[selected] = 2;

                for(int i = 0; i < queryBondListSize; i++)
                {
                    int idx = queryBondList[i];

                    if(query_flags[idx] == 0)
                    {
                        query_flags[idx] = 1;
                        query_parents[idx] = selected;
                    }
                }

                query_order[idx] = selected;
            }
        }


    private:
        bool next_query()
        {
            if(unlikely(core_len >= queryAtomCount))
                return false;

            query_idx = query_order[core_len];
            target_idx = -1;
            target_selector = -1;

            return true;
        }


        bool next_target()
        {
            int query_parent = query_parents[query_idx];

            if(likely(query_parent >= 0))
            {
                int target_parent = core_query[query_parent];
                AtomIdx *targetBondList = molecule_get_bond_list(target, target_parent);
                MolSize targetBondListSize = molecule_get_bond_list_size(target, target_parent);

                for(target_selector++; target_selector < targetBondListSize; target_selector++)
                {
                    int target_id = targetBondList[target_selector];

                    if(!is_core_defined(core_target[target_id]))
                    {
                        target_idx = target_id;
                        return true;
                    }
                }
            }
            else
            {
                for(target_idx++; target_idx < targetAtomCount; target_idx++)
                {
                    if(!is_core_defined(core_target[target_idx]))
                        return true;
                }
            }

            return false;
        }


        bool atom_matches(int queryAtom, int targetAtom)
        {
            int8_t queryAtomNumber = molecule_get_atom_number(query, queryAtom);
            int8_t targetAtomNumber = molecule_get_atom_number(target, targetAtom);

            return queryAtomNumber == targetAtomNumber;
        }


        bool bond_matches(int qIdx1, int qIdx2, int tIdx1, int tIdx2)
        {
            BondIdx queryBond = molecule_get_bond(query, qIdx1, qIdx2);
            BondIdx targetbond = molecule_get_bond(target, tIdx1, tIdx2);

            if(likely(targetbond < 0))
                return false;

            uint8_t queryBondType = molecule_get_bond_type(query, queryBond);
            uint8_t targetbondType = molecule_get_bond_type(target, targetbond);

            return queryBondType == targetbondType;
        }


        bool is_feasible_pair()
        {
            if(likely(!atom_matches(query_idx, target_idx)))
                return false;


            int newQuery = 0;
            int newTarget = 0;

            AtomIdx *queryBondList = molecule_get_bond_list(query, query_idx);
            MolSize queryBondListSize = molecule_get_bond_list_size(query, query_idx);

            for(int i = 0; i < queryBondListSize; i++)
            {
                int other1 = queryBondList[i];

                if(is_core_defined(core_query[other1]))
                {
                    int other2 = core_query[other1];

                    if(!bond_matches(query_idx, other1, target_idx, other2))
                        return false;
                }
                else
                {
                    newQuery++;
                }
            }


            AtomIdx *targetBondList = molecule_get_bond_list(target, target_idx);
            MolSize targetBondListSize = molecule_get_bond_list_size(target, target_idx);

            for(int i = 0; i < targetBondListSize; i++)
            {
                int other2 = targetBondList[i];

                if(!is_core_defined(core_target[other2]))
                    newTarget++;
            }

            return newQuery <= newTarget;
        }


        void undo_add_pair()
        {
            VF2Undo *undo = &undos[--core_len];

            query_idx = query_order[core_len];

            core_query[query_idx] = UNDEFINED_CORE;
            core_target[undo->target_idx] = UNDEFINED_CORE;

            target_selector = undo->target_selector;
            target_idx = undo->target_idx;
        }


        void add_pair()
        {
            VF2Undo *undo = &undos[core_len];

            core_len++;
            core_query[query_idx] = target_idx;
            core_target[target_idx] = query_idx;

            undo->target_selector = target_selector;
            undo->target_idx = target_idx;
        }


        bool is_match_valid()
        {
            std::vector<int> match(core_len);

            for(int i = 0; i < core_len; i++)
                match[i] = core_query[i];

            std::sort(match.begin(), match.end());


            bool included = false;

            for(std::vector<int> &old : matches)
            {
                bool same = true;

                for(int i = 0; i < core_len; i++)
                {
                    if(match[i] != old[i])
                    {
                        same = false;
                        break;
                    }
                }

                if(same)
                {
                    included = true;
                    break;
                }
            }

            if(!included)
                matches.push_back(std::move(match));

            return matches.size() >= limit;
        }


        bool match_core()
        {
            while(true)
            {
                recursion_entry:

                if(unlikely(core_len == query->atomCount))
                {
                    if(is_match_valid())
                        return true;

                    goto recursion_return;
                }

                if(!next_query())
                    goto recursion_return;


                while(next_target())
                {
                    if(unlikely(QueryCancelPending))
                        return false;

                    if(is_feasible_pair())
                    {
                        add_pair();
                        goto recursion_entry;

                        recursion_return:

                        if(core_len == 0)
                            return false;

                        undo_add_pair();
                    }
                }

                goto recursion_return;
            }

            return false;
        }


    public:
        std::vector<std::vector<int>> match(const Molecule *const targetMolecule, int searchLimit)
        {
            target = targetMolecule;
            targetAtomCount = target->atomCount;
            limit = searchLimit;
            core_len = 0;

            std::vector<std::vector<int>> result;

            if(queryAtomCount > targetAtomCount || query->bondCount > target->bondCount)
                return result;


            core_target.resize(targetAtomCount);

            for(int i = 0; i < targetAtomCount; i++)
                core_target[i] = UNDEFINED_CORE;

            for(int i = 0; i < queryAtomCount; i++)
                core_query[i] = UNDEFINED_CORE;

            match_core();

            result.swap(matches);
            return result;
        }
};

#endif /* SUBSTRUCTURE_MATCH_HPP__ */
