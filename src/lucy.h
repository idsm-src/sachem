#ifndef LUCY_H_
#define LUCY_H_

#include <stdlib.h>


typedef struct lucy_Hits Hits;

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


void lucy_init(Lucy *lucy, const char *path);
void lucy_begin(Lucy *lucy);
void lucy_add(Lucy *lucy, int32_t id, const char *fp);
void lucy_delete(Lucy *lucy, int32_t id);
void lucy_commit(Lucy *lucy);
Hits *lucy_search(Lucy *lucy, const char *fp, int max_results);
size_t lucy_get(Lucy *lucy, Hits *hits, int *buffer, size_t size);
void lucy_optimize(Lucy *lucy);

#endif /* LUCY_H_ */
