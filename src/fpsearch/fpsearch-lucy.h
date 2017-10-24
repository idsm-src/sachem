
#ifndef _FPSEARCH_LIB_H
#define _FPSEARCH_LIB_H

#include "../molecule.h"


#ifndef FPSEARCH_API
#define FPSEARCH_API(id) fplucy_##id
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* all parameters for indexing/searching */
typedef struct FPSEARCH_API (params_t)
{
	/* max bonds in subgraph fingerprints */
	unsigned graphSize;

	/* max iterations of ECFP */
	unsigned circSize;

	/* max log-count of indexed feature repetitions */
	unsigned maxLogFeats;

	/* minimum number of overlapping fingerprints on each query atom */
	unsigned searchAtomCoverage;

	/* maximum number of fingerprints in search query */
	unsigned searchMaxFps;

	/* how many index operations to perform before committing to disk */
	unsigned indexerBatch;
}
FPSEARCH_API (params_t);

/* fplucy_initialize
 * returns the initialized `dd` object
 */
void FPSEARCH_API (initialize) (void**ddp, const char *index_dir,
                                const char *fp_ordering_file);

/* fplucy_close
 * safely close the index
 */
void FPSEARCH_API (close) (void *dd);

/* fplucy_params_get
 * load the currently set params to pointed fplucy_params
 */
void FPSEARCH_API (params_get) (void *dd, FPSEARCH_API (params_t) *pp);

/* fplucy_params_set
 * set the parameters from pointed fplucy_params
 */
void FPSEARCH_API (params_set) (void *dd, const FPSEARCH_API (params_t) *pp);

/* fplucy_search
 * returns `ss` search cursor
 */
void* FPSEARCH_API (search) (void *dd,
                             const Molecule *mol,
                             int max_results);

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

/* fplucy_optimize
 * rearranges the index so that searching is faster
 */
void FPSEARCH_API (optimize) (void *dd);

#ifdef __cplusplus
}
#endif

#endif
