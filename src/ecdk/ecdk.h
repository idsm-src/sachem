#ifndef ECDK_H_
#define ECDK_H_

#include <postgres.h>
#include <tcop/tcopprot.h>
#include <stdint.h>
#include "common.h"
#include "java/java.h"


typedef struct
{
    int16_t *counts;
    int16_t *fp;
    uint8_t *molecule;
    bool *restH;

    int fpLength;
} EcdkSubstructureQueryData;


typedef struct
{
    int bitCount;
    ArrayType *fp;
#if USE_COUNT_FINGERPRINT
    ArrayType *counts;
#endif
    bytea *molecule;
    text *error;
} EcdkLoaderData;


void ecdk_java_init(void);
int ecdk_java_parse_substructure_query(EcdkSubstructureQueryData **data, char* query, size_t queryLength, int32_t type, bool implicitHydrogens, bool tautomers);
int ecdk_java_parse_similarity_query(uint64_t **data, char* query, size_t queryLength, int32_t type);
void ecdk_java_parse_data(size_t count, VarChar **molfiles, EcdkLoaderData *data);

#endif /* ECDK_H_ */
