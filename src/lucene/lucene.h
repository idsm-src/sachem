#ifndef LUCENE_H_
#define LUCENE_H_

#include <stdbool.h>
#include <jni.h>
#include "bitset.h"
#include "fingerprints/fingerprint.h"


#define NULL_SUBSEARCH_RESULT       ((LuceneSubsearchResult) { .possition = -1 })


typedef struct
{
    jobject instance;
} Lucene;


typedef struct
{
    jlongArray bitsetArray;
    jlong *bitsetWords;
    BitSet hits;
    size_t possition;
} LuceneSubsearchResult;


typedef struct
{
    jobjectArray result;
    size_t count;
    size_t possition;
} LuceneSimsearchResult;


static inline bool lucene_subsearch_is_open(LuceneSubsearchResult *resultSet)
{
    return resultSet->possition != -1;
}

void lucene_init(Lucene *lucene);
void lucene_terminate(Lucene *lucene);
void lucene_set_folder(Lucene *lucene, const char *path);
LuceneSubsearchResult lucene_subsearch_submit(Lucene *lucene, IntegerFingerprint fp, int32_t maxId);
size_t lucene_subsearch_get(Lucene *lucene, LuceneSubsearchResult *resultSet, int32_t *buffer, size_t size);
void lucene_subsearch_fail(Lucene *lucene, LuceneSubsearchResult *resultSet);
LuceneSimsearchResult lucene_simsearch_submit(Lucene *lucene, IntegerFingerprint fp, int32_t topN, float cutoff);
bool lucene_simsearch_get(Lucene *lucene, LuceneSimsearchResult *result, int32_t *id, float *score);
void lucene_simsearch_fail(Lucene *lucene, LuceneSimsearchResult *result);

#endif /* LUCENE_H_ */
