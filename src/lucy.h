#ifndef LUCY_H_
#define LUCY_H_

#include <stdlib.h>
#include "fingerprints/fingerprint.h"

#define NULL_RESULT_SET     ((LucyResultSet) { .hits = NULL })


typedef struct
{
    struct cfish_String *idF;
    struct cfish_String *fpF;
    struct cfish_String *folder;
    struct lucy_Schema *schema;

    struct lucy_QueryParser *qparser;
    struct lucy_IndexSearcher *searcher;

    struct lucy_Indexer *indexer;
} Lucy;


typedef struct
{
    struct lucy_Hits *hits;
} LucyResultSet;


inline bool lucy_is_open(Lucy *lucy, LucyResultSet *resultSet)
{
    return resultSet->hits != NULL;
}


void lucy_init(Lucy *lucy, const char *path);
void lucy_begin(Lucy *lucy);
void lucy_add(Lucy *lucy, int32_t id, Fingerprint fp);
void lucy_delete(Lucy *lucy, int32_t id);
void lucy_optimize(Lucy *lucy);
void lucy_commit(Lucy *lucy);
void lucy_rollback(Lucy *lucy);
LucyResultSet lucy_search(Lucy *lucy, Fingerprint fp, int max_results);
size_t lucy_get(Lucy *lucy, LucyResultSet *resultSet, int *buffer, size_t size);
void lucy_fail(Lucy *lucy, LucyResultSet *resultSet);

#endif /* LUCY_H_ */
