#ifndef LUCENE_H_
#define LUCENE_H_

#include <stdbool.h>
#include <jni.h>
#include "bitset.h"
#include "fingerprints/fingerprint.h"


#define NULL_RESULT_SET         ((LuceneResultSet) { .possition = -1 })


typedef struct
{
    jobject instance;

    jlongArray bitsetArray;
    jlong *bitsetWords;
    BitSet hits;
} Lucene;


typedef struct
{
    size_t possition;
} LuceneResultSet;


inline bool lucene_is_open(LuceneResultSet *resultSet)
{
    return resultSet->possition != -1;
}


void lucene_init(Lucene *lucene);
void lucene_set_folder(Lucene *lucene, const char *path, int32_t maxId);
void lucene_begin(Lucene *lucene);
void lucene_add(Lucene *lucene, int32_t id, IntegerFingerprint fp);
void lucene_add_index(Lucene *lucene, const char *path);
void lucene_delete(Lucene *lucene, int32_t id);
void lucene_optimize(Lucene *lucene);
void lucene_commit(Lucene *lucene);
void lucene_rollback(Lucene *lucene);
LuceneResultSet lucene_search(Lucene *lucene, IntegerFingerprint fp);
size_t lucene_get(Lucene *lucene, LuceneResultSet *resultSet, int32_t *buffer, size_t size);
void lucene_fail(Lucene *lucene, LuceneResultSet *resultSet);
void lucene_link_directory(const char *oldPath, const char *newPath);
void lucene_delete_directory(const char *path);

#endif /* LUCENE_H_ */
