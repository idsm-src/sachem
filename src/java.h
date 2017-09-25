#ifndef JAVA_H_
#define JAVA_H_

#include <stdbool.h>
#include <stdint.h>


typedef struct
{
    int16_t *counts;
    int16_t *fp;
    uint8_t *atoms;
    uint8_t *bonds;
    bool *restH;

    int fpLength;
    int atomLength;
    int bondLength;
} QueryData;


int java_parse_query(QueryData **data, char* query, size_t queryLength, char *type, bool tautomers);
int java_parse_similarity_query(uint64_t **data, char* query, size_t queryLength, char *type);
void java_module_init(void);
void java_module_finish(void);

#endif /* JAVA_H_ */
