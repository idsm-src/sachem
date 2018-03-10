#ifndef LUCY_H_
#define LUCY_H_

#include <stdbool.h>
#include <stdint.h>
#include "fingerprints/fingerprint.h"

#define USE_ID_TABLE        1
#define NULL_RESULT_SET     ((LucyResultSet) { .possition = -1 })


typedef struct
{
    struct cfish_String *idF;
    struct cfish_String *fpF;
    struct cfish_String *folder;
    struct lucy_Schema *schema;

    struct lucy_QueryParser *qparser;
    struct lucy_IndexSearcher *searcher;
    struct lucy_Collector *collector;
    struct lucy_BitVector *hits;

    struct lucy_Indexer *indexer;

#if USE_ID_TABLE
    int32_t *idTable;
#endif
} Lucy;


typedef struct
{
    size_t possition;
} LucyResultSet;


inline bool lucy_is_open(LucyResultSet *resultSet)
{
    return resultSet->possition != -1;
}


void lucy_init(Lucy *lucy);
void lucy_set_folder(Lucy *lucy, const char *path);
void lucy_begin(Lucy *lucy);
void lucy_add(Lucy *lucy, int32_t id, StringFingerprint fp);
void lucy_add_index(Lucy *lucy, const char *path);
void lucy_delete(Lucy *lucy, int32_t id);
void lucy_optimize(Lucy *lucy);
void lucy_commit(Lucy *lucy);
void lucy_rollback(Lucy *lucy);
LucyResultSet lucy_search(Lucy *lucy, StringFingerprint fp);
size_t lucy_get(Lucy *lucy, LucyResultSet *resultSet, int32_t *buffer, size_t size);
void lucy_fail(Lucy *lucy, LucyResultSet *resultSet);
void lucy_link_directory(const char *oldPath, const char *newPath);
void lucy_delete_directory(const char *path);

#endif /* LUCY_H_ */
