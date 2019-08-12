#ifndef INDEXER_H_
#define INDEXER_H_

#include <stdbool.h>
#include <jni.h>
#include "bitset.h"
#include "fingerprints/fingerprint.h"


typedef struct
{
    jobject instance;
} LuceneIndexer;


void lucene_indexer_init(LuceneIndexer *lucene);
void lucene_indexer_terminate(LuceneIndexer *lucene);
void lucene_indexer_begin(LuceneIndexer *lucene, const char *path);
void lucene_indexer_add(LuceneIndexer *lucene, int32_t id, IntegerFingerprint subfp, IntegerFingerprint simfp);
void lucene_indexer_add_index(LuceneIndexer *lucene, const char *path);
void lucene_indexer_delete(LuceneIndexer *lucene, int32_t id);
void lucene_indexer_optimize(LuceneIndexer *lucene);
void lucene_indexer_commit(LuceneIndexer *lucene);
void lucene_indexer_rollback(LuceneIndexer *lucene);
void lucene_indexer_link_directory(const char *oldPath, const char *newPath);
void lucene_indexer_delete_directory(const char *path);

#endif /* INDEXER_H_ */
