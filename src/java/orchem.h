#ifndef JAVA_ORCHEM_H_
#define JAVA_ORCHEM_H_

#include <postgres.h>
#include <tcop/tcopprot.h>
#include <stdint.h>
#include "java.h"


typedef struct
{
    int16_t *counts;
    int16_t *fp;
    uint8_t *molecule;
    bool *restH;

    int fpLength;
} OrchemSubstructureQueryData;


typedef struct
{
    int bitCount;
    ArrayType *fp;
    ArrayType *counts;
    bytea *molecule;
    text *error;
} OrchemLoaderData;


void java_orchem_init(void);
int java_orchem_parse_substructure_query(OrchemSubstructureQueryData **data, char* query, size_t queryLength, int32_t type, bool implicitHydrogens, bool tautomers);
int java_orchem_parse_similarity_query(uint64_t **data, char* query, size_t queryLength, int32_t type);
void java_orchem_parse_data(size_t count, VarChar **molfiles, OrchemLoaderData *data);

#endif /* JAVA_ORCHEM_H_ */
