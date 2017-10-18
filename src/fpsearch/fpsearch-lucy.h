
#ifndef _FPSEARCH_LIB_H
#define _FPSEARCH_LIB_H

#define restrict //lol
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)
#define MOLECULE_H_NO_POSTGRES
#include "../molecule.h"
#undef restrict

#define FPSEARCH_API(id) fplucy_##id

#ifdef __cplusplus
extern "C" {
#endif

/* fplucy_initialize
 * returns the initialized `dd` object
 */
void* FPSEARCH_API (initialize) (const char *index_dir,
                                 const char *fp_ordering_file);

/* fplucy_close
 * safely close the index
 */
void FPSEARCH_API (close) (void *dd);

/* fplucy_search
 * returns `ss` search cursor
 */
void* FPSEARCH_API (search) (void *dd,
                             const Molecule *mol,
                             int max_results = 1000000);

/* fplucy_search_fillbuf
 * output `n` matches of `ss` to pre-allocated array `results`,
 * returning the amount of results retrieved
 */
int FPSEARCH_API (search_fillbuf) (void *dd,
                                   void *ss,
                                   int *results,
                                   int n);

/* fplucy_search_finish
 * free the `ss` search cursor
 */
void FPSEARCH_API (search_finish) (void *dd,
                                   void *ss);

/* fplucy_add_mol
 * adds molecule with `guid` identifier to index `dd`
 */
void FPSEARCH_API (add_mol) (void *dd,
                             int guid,
                             const Molecule*mol);

/* fplucy_remove_mol
 * removes the molecule identified by `guid` from the index `dd`
 */
void FPSEARCH_API (remove_mol) (void *dd,
                                int guid);

#ifdef __cplusplus
}
#endif

#endif
