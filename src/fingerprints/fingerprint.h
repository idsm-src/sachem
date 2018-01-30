#ifndef FINGERPRINT_H__
#define FINGERPRINT_H__

#include <stdlib.h>


typedef struct Molecule Molecule;

typedef struct
{
    size_t size;
    char *data;
} Fingerprint;


Fingerprint fingerprint_get(const Molecule *molecule, void *(*alloc)(size_t));
Fingerprint fingerprint_get_query(const Molecule *molecule, void *(*alloc)(size_t));


inline bool fingerprint_is_valid(Fingerprint fingerprint)
{
    return fingerprint.size != (size_t) -1;
}

#endif
