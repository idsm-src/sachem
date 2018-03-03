#ifndef FINGERPRINT_H__
#define FINGERPRINT_H__

#include <stdbool.h>
#include <stdlib.h>


#define GRAPH_SIZE              7
#define RING_SIZE               3
#define MAX_FEAT_LOGCOUNT       5


typedef struct Molecule Molecule;


typedef struct
{
    size_t size;
    char *data;
} StringFingerprint;


typedef struct
{
    size_t size;
    int32_t *data;
} IntegerFingerprint;


StringFingerprint string_fingerprint_get(const Molecule *molecule, void *(*alloc)(size_t));
StringFingerprint string_fingerprint_get_query(const Molecule *molecule, void *(*alloc)(size_t));

IntegerFingerprint integer_fingerprint_get(const Molecule *molecule, void *(*alloc)(size_t));
IntegerFingerprint integer_fingerprint_get_query(const Molecule *molecule, void *(*alloc)(size_t));


inline bool string_fingerprint_is_valid(StringFingerprint fingerprint)
{
    return fingerprint.size != (size_t) -1;
}


inline bool integer_fingerprint_is_valid(IntegerFingerprint fingerprint)
{
    return fingerprint.size != (size_t) -1;
}

#endif
