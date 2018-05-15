#ifndef JAVA_PARSE_H_
#define JAVA_PARSE_H_

#include <postgres.h>
#include <tcop/tcopprot.h>
#include <stdint.h>
#include "java.h"


typedef struct
{
    uint8_t *molecule;
    bool *restH;
} SubstructureQueryData;


typedef struct
{
    uint8_t *molecule;
} SimilarityQueryData;


typedef struct
{
    bytea *molecule;
    text *error;
} LoaderData;


void java_parse_init(void);
int java_parse_substructure_query(SubstructureQueryData **data, char* query, size_t queryLength, int32_t type, bool implicitHydrogens, bool tautomers);
void java_parse_similarity_query(SimilarityQueryData *data, char* query, size_t queryLength, int32_t type);
void java_parse_data(size_t count, VarChar **molfiles, LoaderData *data);

#endif /* JAVA_PARSE_H_ */
