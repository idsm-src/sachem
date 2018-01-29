#ifndef FINGERPRINT_H__
#define FINGERPRINT_H__

#include <stdlib.h>


typedef struct Molecule Molecule;


char *fingerprint_get(const Molecule *molecule, void *(*alloc)(size_t));
char *fingerprint_get_query(const Molecule *molecule, void *(*alloc)(size_t));

#endif
