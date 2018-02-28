#ifndef LUCENE_H_
#define LUCENE_H_

#include <stdbool.h>
#include <jni.h>
#include "fingerprints/fingerprint.h"

#define NULL_RESULT_SET     NULL


typedef jobject Lucene;
typedef jobject LuceneResultSet;


inline bool lucene_is_open(Lucene *lucene, LuceneResultSet *resultSet)
{
    return *resultSet != NULL;
}


void lucene_init(Lucene *lucene);
void lucene_set_folder(Lucene *lucene, const char *path);
void lucene_begin(Lucene *lucene);
void lucene_add(Lucene *lucene, int32_t id, IntegerFingerprint fp);
void lucene_delete(Lucene *lucene, int32_t id);
void lucene_optimize(Lucene *lucene);
void lucene_commit(Lucene *lucene);
void lucene_rollback(Lucene *lucene);
LuceneResultSet lucene_search(Lucene *lucene, IntegerFingerprint fp, int max_results);
size_t lucene_get(Lucene *lucene, LuceneResultSet *resultSet, int *buffer, size_t size);
void lucene_fail(Lucene *lucene, LuceneResultSet *resultSet);
void lucene_link_directory(const char *oldPath, const char *newPath);
void lucene_delete_directory(const char *path);

#endif /* LUCENE_H_ */