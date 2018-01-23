#ifndef JAVA_H_
#define JAVA_H_

#include <postgres.h>
#include <utils/array.h>
#include <stdbool.h>
#include <stdint.h>


typedef struct
{
    uint8_t *molecule;
    bool *restH;
} SubstructureQueryData;


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


int java_parse_substructure_query(SubstructureQueryData **data, char* query, size_t queryLength, int32_t type, bool implicitHydrogens, bool tautomers);
int java_parse_orchem_substructure_query(OrchemSubstructureQueryData **data, char* query, size_t queryLength, int32_t type, bool implicitHydrogens, bool tautomers);
int java_parse_orchem_similarity_query(uint64_t **data, char* query, size_t queryLength, int32_t type);
void java_parse_orchem_data(size_t count, VarChar **molfiles, OrchemLoaderData *data);
void java_module_init(void);
void java_module_finish(void);

#endif /* JAVA_H_ */
