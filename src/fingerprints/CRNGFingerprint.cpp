#include "CRNGFingerprint.hpp"
#include "SubstructureMatch.hpp"

extern "C"
{
#include <postgres.h>
#include <utils/memutils.h>
#include "sachem.h"
#include "molecule.h"
#include "isomorphism.h"
}


#define PATTERN_COUNT 345


static bool initialized = false;
static Molecule patternMolecule[PATTERN_COUNT];

static uint8_t *patterns[PATTERN_COUNT] = {
    (uint8_t []) {0,0,0,10,0,0,0,10,0,2,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,2,6,0,7,1,7,0,8,1,8,0,9,1,0,0,9,1,48,0,3,48,5,3},
    (uint8_t []) {0,0,0,12,0,0,0,12,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,8,0,9,1,9,0,10,1,10,0,11,1,0,0,11,1,48,0,3},
    (uint8_t []) {0,1,0,11,0,0,0,12,0,0,8,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,0,1,0,0,7,1,7,0,8,1,8,0,9,1,9,0,10,1,10,0,11,1,1,0,11,1},
    (uint8_t []) {0,3,0,16,0,0,0,19,0,0,8,8,8,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,8,0,0,1,0,0,9,1,9,0,10,1,10,0,11,1,11,0,1,1,1,0,12,1,12,0,13,1,13,0,2,1,2,0,14,1,14,0,15,1,15,0,16,1,16,0,17,1,17,0,18,1,3,0,18,1},
    (uint8_t []) {0,1,0,6,0,0,0,7,0,1,8,1,0,2,2,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,1,0,6,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,2,7,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,1,0,0,1,48,0,3},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,7,15,8,3,0,4,1,4,0,0,1,0,0,1,1,1,0,2,1,2,0,5,1,3,0,5,1},
    (uint8_t []) {0,3,0,9,0,0,0,12,0,0,8,8,8,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,3,0,11,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,8,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1,0,1,1},
    (uint8_t []) {0,8,0,11,0,0,0,19,0,2,15,27,8,8,7,7,7,7,8,0,4,11,4,0,9,1,9,0,10,1,10,0,11,1,11,0,2,1,2,0,0,1,0,0,3,1,3,0,12,1,12,0,13,1,13,0,5,1,5,0,14,1,14,0,15,1,15,0,16,1,16,0,17,1,17,0,18,1,18,0,6,1,6,0,1,1,1,0,7,1,8,0,7,11,0,1,253,0,7,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,16,2,0,0,11,0,0,3,1,3,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,1,7,7,7,3,0,0,2,0,0,4,1,4,0,1,1,1,0,5,1,5,0,2,1,3,0,2,1,48,0,3},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,7,7,7,3,0,0,1,0,0,4,1,4,0,1,1,1,0,2,1,3,0,2,1},
    (uint8_t []) {0,6,0,0,0,0,0,6,0,0,8,74,8,74,8,74,0,0,1,1,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,0,0,5,1},
    (uint8_t []) {0,0,0,4,0,0,0,4,0,0,0,0,1,11,1,0,2,11,2,0,3,11,0,0,3,11},
    (uint8_t []) {0,0,0,9,0,0,0,9,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,0,0,8,1,48,0,3},
    (uint8_t []) {0,0,0,11,0,0,0,11,0,2,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,8,2,8,0,9,1,9,0,10,1,0,0,10,1,48,0,3,48,7,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,1,2,0,3,11,3,0,4,1,4,0,5,1,5,0,0,1,1,0,0,2,48,5,3},
    (uint8_t []) {0,2,0,5,0,0,0,7,0,1,7,7,2,0,3,1,3,0,4,11,4,0,0,1,0,0,5,11,5,0,6,1,6,0,1,1,2,0,1,2,48,6,3},
    (uint8_t []) {0,4,0,12,0,0,0,16,0,0,7,7,7,7,4,0,5,1,5,0,6,1,6,0,0,11,0,0,7,11,7,0,8,1,8,0,9,1,9,0,1,11,1,0,10,11,10,0,11,1,11,0,12,1,12,0,2,11,2,0,13,11,13,0,14,1,14,0,15,1,15,0,3,11,4,0,3,11},
    (uint8_t []) {0,1,0,8,0,0,0,9,0,1,8,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,1,0,8,1,48,0,3},
    (uint8_t []) {0,4,0,24,0,0,0,28,0,0,8,7,8,7,4,0,5,11,5,0,6,11,6,0,0,1,0,0,7,1,7,0,8,11,8,0,9,11,9,0,10,11,10,0,11,1,11,0,12,1,12,0,13,1,13,0,1,1,1,0,14,1,14,0,15,1,15,0,16,1,16,0,17,11,17,0,18,11,18,0,19,11,19,0,2,1,2,0,20,1,20,0,21,11,21,0,22,11,22,0,23,1,23,0,24,1,24,0,3,1,3,0,25,1,25,0,26,1,26,0,27,1,4,0,27,1},
    (uint8_t []) {0,7,0,13,0,0,0,20,0,0,16,16,7,7,7,7,7,7,0,8,1,8,0,2,1,2,0,9,1,9,0,10,1,10,0,3,1,3,0,11,1,11,0,12,1,12,0,4,1,4,0,13,1,13,0,14,1,14,0,5,1,5,0,15,1,15,0,16,1,16,0,6,1,6,0,17,1,17,0,18,1,18,0,0,1,0,0,1,1,1,0,19,1,7,0,19,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,1,3,0,0,11,0,0,1,11,1,0,4,11,2,0,4,11},
    (uint8_t []) {0,1,0,6,0,0,0,7,0,0,8,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,11,4,0,5,1,5,0,6,1,1,0,6,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,2,7,7,2,0,0,1,0,0,3,1,3,0,4,1,4,0,5,1,5,0,1,1,2,0,1,2,0,1,1,48,5,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,2,7,7,2,0,0,1,0,0,3,1,3,0,4,1,4,0,1,1,2,0,1,2,0,1,1,48,4,3},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,16,7,7,3,0,1,1,1,0,4,1,4,0,0,1,0,0,5,1,5,0,2,1,3,0,2,1},
    (uint8_t []) {0,8,0,0,0,0,0,8,0,0,8,74,8,74,8,74,8,74,0,0,1,1,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,0,0,7,1},
    (uint8_t []) {0,0,0,8,0,0,0,8,0,0,0,0,1,11,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,0,0,7,1},
    (uint8_t []) {0,0,0,13,0,0,0,13,0,3,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,2,5,0,6,1,6,0,7,1,7,0,8,1,8,0,9,1,9,0,10,2,10,0,11,1,11,0,12,1,0,0,12,1,48,0,3,48,4,3,48,9,3},
    (uint8_t []) {0,3,0,26,0,0,0,29,0,6,7,8,8,3,0,4,2,4,0,5,1,5,0,6,2,6,0,7,1,7,0,8,1,8,0,9,1,9,0,1,1,1,0,10,1,10,0,11,1,11,0,12,1,12,0,0,1,0,0,13,1,13,0,14,1,14,0,2,1,2,0,15,1,15,0,16,1,16,0,17,1,17,0,18,1,18,0,19,1,19,0,20,2,20,0,21,1,21,0,22,1,22,0,23,1,23,0,24,1,24,0,25,1,25,0,26,1,26,0,27,1,27,0,28,2,3,0,28,1,48,0,3,48,1,3,48,2,3,48,19,3,48,27,3,48,28,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,11,5,0,0,1,1,0,0,2,48,5,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,3,8,1,0,2,1,2,0,3,2,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,2,48,0,3,48,1,3,48,5,3},
    (uint8_t []) {0,14,0,26,0,0,0,40,0,0,8,8,8,8,8,8,8,8,8,8,8,8,8,8,14,0,15,1,15,0,16,1,16,0,0,1,0,0,17,1,17,0,1,1,1,0,18,1,18,0,19,1,19,0,2,1,2,0,20,1,20,0,3,1,3,0,21,1,21,0,22,1,22,0,4,1,4,0,23,1,23,0,5,1,5,0,24,1,24,0,25,1,25,0,6,1,6,0,26,1,26,0,27,1,27,0,28,1,28,0,29,1,29,0,7,1,7,0,30,1,30,0,8,1,8,0,31,1,31,0,32,1,32,0,9,1,9,0,33,1,33,0,10,1,10,0,34,1,34,0,35,1,35,0,11,1,11,0,36,1,36,0,12,1,12,0,37,1,37,0,38,1,38,0,13,1,13,0,39,1,14,0,39,1},
    (uint8_t []) {0,4,0,11,0,0,0,15,0,11,7,7,7,7,4,0,5,2,5,0,0,1,0,0,6,2,6,0,7,1,7,0,8,2,8,0,1,1,1,0,9,1,9,0,10,1,10,0,2,1,2,0,11,2,11,0,12,1,12,0,13,2,13,0,3,1,3,0,14,2,4,0,14,1,48,0,3,48,1,3,48,2,3,48,3,3,48,4,3,48,9,3,48,10,3,48,11,3,48,12,3,48,13,3,48,14,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,1,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,7,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,2,48,4,3},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,5,27,7,7,3,0,4,2,4,0,1,1,1,0,0,1,0,0,2,1,2,0,5,2,3,0,5,1,0,0,253,0,2,1,48,0,3,48,4,3,48,5,3},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,2,27,7,7,3,0,4,1,4,0,1,1,1,0,0,1,0,0,2,1,3,0,2,1,0,0,253,0,1,1},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,7,7,7,3,0,4,11,4,0,0,11,0,0,1,11,1,0,5,11,5,0,2,11,3,0,2,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,16,7,2,0,3,11,3,0,0,1,0,0,1,1,1,0,4,1,2,0,4,1},
    (uint8_t []) {0,4,0,1,0,0,0,5,0,0,7,7,7,7,4,0,0,11,0,0,1,11,1,0,2,11,2,0,3,11,4,0,3,11},
    (uint8_t []) {0,4,0,0,0,0,0,4,0,1,16,26,16,26,0,0,1,1,1,0,2,1,2,0,3,1,0,0,3,1,0,3,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,11,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,11,1,0,5,1,0,0,1},
    (uint8_t []) {0,3,0,11,0,0,0,14,0,0,8,7,7,3,0,4,11,4,0,5,11,5,0,0,1,0,0,6,1,6,0,7,11,7,0,8,11,8,0,9,1,9,0,1,1,1,0,10,1,10,0,11,1,11,0,2,1,2,0,12,1,12,0,13,1,3,0,13,1},
    (uint8_t []) {0,3,0,12,0,0,0,15,0,0,8,8,8,3,0,4,11,4,0,5,1,5,0,0,1,0,0,6,1,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,11,0,12,1,12,0,13,1,13,0,14,1,3,0,14,1},
    (uint8_t []) {0,2,0,5,0,0,0,7,0,0,7,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,1,0,6,1,2,0,6,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,3,1,3,0,0,11,0,0,4,11,4,0,1,11,1,0,5,11,2,0,5,11,0,1,255},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,2,7,7,2,0,3,2,3,0,0,1,0,0,4,1,4,0,1,2,1,0,5,1,2,0,5,1,48,0,3,48,3,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,8,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,1,0,5,11,2,0,5,11},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,3,27,7,7,3,0,4,2,4,0,1,1,1,0,0,1,0,0,2,1,2,0,5,11,3,0,5,1,0,0,254,0,1,1,48,0,3},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,7,7,7,3,0,4,11,4,0,0,11,0,0,1,11,1,0,2,11,2,0,5,11,3,0,5,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,16,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,11,4,0,5,1,1,0,5,1},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,16,7,7,3,0,4,11,4,0,0,1,0,0,1,1,1,0,5,1,5,0,2,1,3,0,2,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,0,2,0,0,3,1,3,0,4,2,4,0,5,1,5,0,1,1,2,0,1,11,48,2,3},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,11,4,0,5,11,0,0,5,1,48,0,3},
    (uint8_t []) {0,0,0,7,0,0,0,7,0,0,0,0,1,1,1,0,2,11,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,0,0,6,11},
    (uint8_t []) {0,0,0,10,0,0,0,10,0,0,0,0,1,1,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,8,0,9,1,0,0,9,1},
    (uint8_t []) {0,0,0,14,0,0,0,14,0,2,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,8,0,9,1,9,0,10,1,10,0,11,2,11,0,12,1,12,0,13,1,0,0,13,1,48,0,3,48,10,3},
    (uint8_t []) {0,1,0,7,0,0,0,8,0,0,7,1,0,2,11,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,1,0,7,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,4,11,4,0,0,11,0,0,1,11,1,0,5,11,2,0,5,11},
    (uint8_t []) {0,4,0,8,0,0,0,12,0,0,7,7,7,7,4,0,5,1,5,0,0,1,0,0,6,1,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,2,1,2,0,10,1,10,0,11,1,11,0,3,1,4,0,3,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,3,7,7,2,0,3,1,3,0,0,1,0,0,4,2,4,0,1,1,2,0,1,2,48,2,3,48,3,3,48,4,3},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,1,7,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,1,1,0,4,1,0,0,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11,0,0,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,16,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,1,1,0,4,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,16,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,5,0,1,11,2,0,1,11,0,0,1},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,3,5,7,7,0,0,1,1,1,0,3,11,3,0,4,1,4,0,5,2,5,0,2,1,0,0,2,1,0,0,255,0,2,1,48,3,3},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,0,0,0,1,1,1,0,2,11,2,0,3,11,3,0,4,11,0,0,4,11},
    (uint8_t []) {0,0,0,11,0,0,0,11,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,8,0,9,1,9,0,10,1,0,0,10,1,48,0,3},
    (uint8_t []) {0,1,0,8,0,0,0,9,0,0,7,1,0,2,11,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,1,0,8,1},
    (uint8_t []) {0,2,0,9,0,0,0,11,0,0,8,8,2,0,3,1,3,0,4,11,4,0,5,1,5,0,0,1,0,0,6,1,6,0,7,1,7,0,8,1,8,0,1,1,1,0,9,1,9,0,10,1,2,0,10,11},
    (uint8_t []) {0,4,0,6,0,0,0,10,0,3,12,7,7,7,4,0,5,2,5,0,6,1,6,0,1,2,1,0,7,1,7,0,8,2,8,0,9,1,9,0,2,11,2,0,0,1,0,0,3,1,4,0,3,11,48,2,3,48,3,3,48,4,3},
    (uint8_t []) {0,5,0,11,0,0,0,16,0,0,8,7,7,7,7,5,0,6,1,6,0,1,1,1,0,7,1,7,0,8,1,8,0,2,1,2,0,9,1,9,0,10,1,10,0,3,1,3,0,11,1,11,0,12,1,12,0,4,1,4,0,13,1,13,0,14,1,14,0,0,1,0,0,15,1,5,0,15,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,3,64,7,7,3,0,4,1,4,0,1,1,1,0,0,1,0,0,2,1,3,0,2,1,0,0,255,0,1,1,0,2,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,5,0,1,11,2,0,1,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,16,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,16,2,0,0,11,0,0,3,1,3,0,4,1,4,0,1,1,2,0,1,1},
    (uint8_t []) {0,4,0,1,0,0,0,5,0,1,7,7,7,7,4,0,0,11,0,0,1,11,1,0,2,11,2,0,3,11,4,0,3,11,0,2,1},
    (uint8_t []) {0,8,0,0,0,0,0,8,0,0,8,14,8,14,8,14,8,14,0,0,1,1,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,0,0,7,1},
    (uint8_t []) {0,0,0,8,0,0,0,8,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,0,0,7,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,2,7,1,0,2,11,2,0,3,1,3,0,4,2,4,0,5,1,5,0,0,1,1,0,0,2,0,0,1,48,2,3},
    (uint8_t []) {0,4,0,15,0,0,0,19,0,1,8,7,7,7,4,0,5,2,5,0,6,1,6,0,7,1,7,0,8,1,8,0,0,1,0,0,9,1,9,0,10,1,10,0,11,1,11,0,1,1,1,0,12,1,12,0,13,1,13,0,2,1,2,0,14,1,14,0,15,1,15,0,3,1,3,0,16,1,16,0,17,1,17,0,18,1,4,0,18,1,48,0,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11,0,0,255},
    (uint8_t []) {0,3,0,4,0,0,0,7,0,0,8,7,7,3,0,4,11,4,0,1,11,1,0,5,1,5,0,0,1,0,0,6,1,6,0,2,1,3,0,2,11},
    (uint8_t []) {0,4,0,11,0,0,0,15,0,0,8,8,7,7,4,0,5,11,5,0,0,1,0,0,6,1,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,11,0,12,1,12,0,13,1,13,0,3,1,3,0,14,1,4,0,14,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,16,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,1,4,0,5,1,1,0,5,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,16,7,7,3,0,1,11,1,0,4,11,4,0,0,11,0,0,2,11,3,0,2,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,11,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,11,1,0,5,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,1,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,11},
    (uint8_t []) {0,7,0,16,0,0,0,23,0,0,7,7,7,7,7,7,7,7,0,8,1,8,0,9,1,9,0,0,1,0,0,10,1,10,0,11,1,11,0,1,1,1,0,12,1,12,0,13,1,13,0,2,1,2,0,14,1,14,0,15,1,15,0,3,1,3,0,16,1,16,0,17,1,17,0,4,1,4,0,18,1,18,0,19,1,19,0,5,1,5,0,20,1,20,0,21,1,21,0,6,1,6,0,22,1,7,0,22,1},
    (uint8_t []) {0,2,0,5,0,0,0,7,0,0,7,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,11,5,0,1,1,1,0,6,1,2,0,6,1},
    (uint8_t []) {0,2,0,5,0,0,0,7,0,0,8,8,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,11,5,0,1,1,1,0,6,1,2,0,6,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,16,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,0,0,0,1,11,1,0,2,1,2,0,3,1,3,0,4,11,4,0,5,11,0,0,5,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,5,0,1,11,2,0,1,11,0,0,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,0,11,1,0,0,11},
    (uint8_t []) {0,1,0,7,0,0,0,8,0,0,8,1,0,2,1,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,6,0,7,1,1,0,7,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,11,2,0,3,11,3,0,0,1,0,0,4,1,1,0,4,1},
    (uint8_t []) {0,1,0,6,0,0,0,7,0,0,7,1,0,2,11,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,6,1,1,0,6,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,2,2,0,0,1,0,0,3,1,3,0,4,11,4,0,5,1,1,0,5,1,48,0,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,8,2,0,3,1,3,0,0,1,0,0,4,1,4,0,1,1,1,0,5,1,2,0,5,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,8,2,0,3,1,3,0,0,1,0,0,1,1,1,0,4,1,2,0,4,1},
    (uint8_t []) {0,2,0,16,0,0,0,18,0,0,8,8,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,11,5,0,6,11,6,0,7,1,7,0,8,1,8,0,9,1,9,0,10,11,10,0,11,11,11,0,1,1,1,0,12,1,12,0,13,11,13,0,14,11,14,0,15,11,15,0,16,1,16,0,17,1,2,0,17,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,8,8,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,11,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,16,7,2,0,3,2,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,1,48,0,3},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,1,16,7,7,3,0,4,1,4,0,0,1,0,0,5,1,5,0,1,11,1,0,2,1,3,0,2,2,48,5,3},
    (uint8_t []) {0,1,0,2,0,0,0,3,0,0,40,1,0,2,1,2,0,0,1,1,0,0,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,0,2,0,0,3,1,3,0,4,11,4,0,5,1,5,0,1,1,2,0,1,1,48,0,3},
    (uint8_t []) {0,0,0,11,0,0,0,11,0,2,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,2,7,0,8,1,8,0,9,1,9,0,10,1,0,0,10,1,48,0,3,48,6,3},
    (uint8_t []) {0,1,0,8,0,0,0,9,0,0,8,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,0,1,0,0,6,1,6,0,7,1,7,0,8,1,1,0,8,1},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,1,0,0,1,11,1,0,2,1,2,0,3,1,3,0,4,1,0,0,4,1,0,3,255},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,11,2,0,3,1,3,0,0,1,0,0,4,11,1,0,4,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,3,7,1,0,2,2,2,0,0,1,0,0,3,2,3,0,4,1,1,0,4,1,48,0,3,48,1,3,48,2,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1,0,0,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,8,8,2,0,3,1,3,0,0,1,0,0,1,1,1,0,4,1,2,0,4,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,1,7,7,7,3,0,0,2,0,0,4,1,4,0,1,1,1,0,2,1,3,0,2,1,48,0,3},
    (uint8_t []) {0,4,0,1,0,0,0,5,0,2,7,26,8,8,4,0,0,1,0,0,2,1,2,0,1,1,1,0,3,1,4,0,3,2,0,1,253,0,3,1},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,0,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,0,0,4,11},
    (uint8_t []) {0,2,0,10,0,0,0,12,0,0,8,8,2,0,3,1,3,0,4,1,4,0,5,1,5,0,0,1,0,0,6,1,6,0,7,1,7,0,8,1,8,0,9,1,9,0,1,1,1,0,10,1,10,0,11,1,2,0,11,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,8,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,5,0,1,11,2,0,1,11,0,0,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,16,7,7,3,0,4,11,4,0,0,11,0,0,1,11,1,0,2,11,3,0,2,11},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,7,7,7,3,0,0,1,0,0,4,1,4,0,1,1,1,0,5,1,5,0,2,1,3,0,2,1},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,0,0,0,1,11,1,0,2,1,2,0,3,11,3,0,4,11,0,0,4,1},
    (uint8_t []) {0,2,0,22,0,0,0,24,0,4,7,8,2,0,3,1,3,0,4,2,4,0,5,1,5,0,0,1,0,0,6,1,6,0,7,11,7,0,8,11,8,0,9,11,9,0,10,11,10,0,11,1,11,0,12,1,12,0,1,1,1,0,13,1,13,0,14,2,14,0,15,1,15,0,16,1,16,0,17,1,17,0,18,1,18,0,19,1,19,0,20,1,20,0,21,1,21,0,22,1,22,0,23,1,2,0,23,2,48,0,3,48,1,3,48,13,3,48,23,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,11,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,16,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,1,42,16,16,3,0,4,2,4,0,1,1,1,0,0,1,0,0,2,1,3,0,2,1,48,0,3},
    (uint8_t []) {0,3,0,11,0,0,0,14,0,0,8,8,8,3,0,4,11,4,0,5,1,5,0,0,1,0,0,6,1,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,11,0,12,1,12,0,13,1,3,0,13,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,2,27,7,7,3,0,4,1,4,0,1,1,1,0,0,1,0,0,2,1,3,0,2,1,0,0,254,0,1,1},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,6,27,7,7,3,0,4,2,4,0,1,1,1,0,0,1,0,0,2,1,2,0,5,2,3,0,5,1,0,0,253,0,1,1,0,2,1,48,0,3,48,4,3,48,5,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,1,3,0,0,1,0,0,1,1,1,0,4,1,2,0,4,1},
    (uint8_t []) {0,4,0,1,0,0,0,5,0,0,7,7,7,7,4,0,0,11,0,0,1,11,1,0,2,11,2,0,3,11,4,0,3,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,8,1,0,2,1,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,1},
    (uint8_t []) {0,4,0,10,0,0,0,14,0,0,8,8,7,7,4,0,5,11,5,0,0,1,0,0,6,1,6,0,7,1,7,0,8,1,8,0,1,1,1,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,11,0,12,1,12,0,3,1,3,0,13,1,4,0,13,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,1,2,0,0,11,0,0,3,11,3,0,4,11,4,0,5,11,1,0,5,11},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,3,7,7,2,0,3,2,3,0,0,1,0,0,4,2,4,0,1,1,1,0,5,1,2,0,5,1,48,0,3,48,1,3,48,2,3},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,6,27,7,7,3,0,4,2,4,0,1,1,1,0,0,1,0,0,2,1,2,0,5,2,3,0,5,1,0,0,254,0,1,1,0,2,1,48,0,3,48,4,3,48,5,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11},
    (uint8_t []) {0,2,0,16,0,0,0,18,0,0,8,8,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,11,5,0,6,11,6,0,7,1,7,0,8,1,8,0,9,1,9,0,10,11,10,0,11,11,11,0,12,11,12,0,1,1,1,0,13,1,13,0,14,11,14,0,15,11,15,0,16,1,16,0,17,1,2,0,17,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,8,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11,0,1,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,16,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,15,0,25,0,0,0,40,0,0,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,15,0,16,1,16,0,17,1,17,0,0,1,0,0,18,1,18,0,1,1,1,0,19,1,19,0,20,1,20,0,2,1,2,0,21,1,21,0,3,1,3,0,22,1,22,0,23,1,23,0,4,1,4,0,24,1,24,0,5,1,5,0,25,1,25,0,26,1,26,0,6,1,6,0,27,1,27,0,7,1,7,0,28,1,28,0,29,1,29,0,8,1,8,0,30,1,30,0,9,1,9,0,31,1,31,0,32,1,32,0,10,1,10,0,33,1,33,0,11,1,11,0,34,1,34,0,35,1,35,0,12,1,12,0,36,1,36,0,13,1,13,0,37,1,37,0,38,1,38,0,14,1,14,0,39,1,15,0,39,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,4,7,1,0,2,2,2,0,0,1,0,0,3,2,3,0,4,1,1,0,4,1,0,0,1,48,0,3,48,1,3,48,2,3},
    (uint8_t []) {0,2,0,5,0,0,0,7,0,0,8,8,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,5,0,6,1,6,0,1,1,2,0,1,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,8,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,1},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,11,0,0,4,1,48,0,3},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,3,0,0,1,1,1,0,2,2,2,0,3,1,3,0,4,1,0,0,4,2,48,0,3,48,1,3,48,4,3},
    (uint8_t []) {0,1,0,6,0,0,0,7,0,0,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,0,1,1,0,0,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,11,3,0,0,11,0,0,1,11,1,0,4,11,2,0,4,11},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,7,16,16,3,0,4,1,4,0,1,1,1,0,2,1,2,0,5,1,5,0,0,1,3,0,0,1},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,0,0,0,1,11,1,0,2,11,2,0,3,1,3,0,4,1,0,0,4,1},
    (uint8_t []) {0,0,0,7,0,0,0,7,0,0,0,0,1,11,1,0,2,11,2,0,3,11,3,0,4,11,4,0,5,11,5,0,6,11,0,0,6,11},
    (uint8_t []) {0,1,0,7,0,0,0,8,0,0,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,6,0,7,1,1,0,7,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,8,1,0,2,2,2,0,3,1,3,0,4,11,4,0,0,1,0,0,5,1,1,0,5,11},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,11,2,0,0,2,0,0,3,1,3,0,4,1,1,0,4,2},
    (uint8_t []) {0,1,0,6,0,0,0,7,0,0,7,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,11,4,0,5,1,5,0,6,1,1,0,6,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,1,3,0,0,11,0,0,4,11,4,0,5,11,5,0,1,11,2,0,1,11},
    (uint8_t []) {0,3,0,8,0,0,0,11,0,0,8,7,7,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,0,1,0,0,10,1,10,0,2,1,3,0,2,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,3,7,1,0,2,1,2,0,3,2,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,2,48,0,3,48,1,3,48,5,3},
    (uint8_t []) {0,2,0,5,0,0,0,7,0,1,7,7,2,0,3,1,3,0,4,11,4,0,0,1,0,0,5,1,5,0,6,1,6,0,1,1,2,0,1,2,48,6,3},
    (uint8_t []) {0,4,0,10,0,0,0,14,0,0,8,8,7,7,4,0,5,11,5,0,0,1,0,0,6,1,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,11,0,12,1,12,0,3,1,3,0,13,1,4,0,13,1},
    (uint8_t []) {0,1,0,2,0,0,0,3,0,0,92,1,0,2,1,2,0,0,1,1,0,0,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,11,2,0,3,1,3,0,4,1,4,0,5,1,5,0,0,1,1,0,0,2,0,0,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,8,7,2,0,3,1,3,0,0,11,0,0,4,11,4,0,5,11,5,0,1,11,2,0,1,11},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,15,8,8,3,0,4,1,4,0,1,1,1,0,0,1,0,0,2,1,3,0,2,1},
    (uint8_t []) {0,0,0,8,0,0,0,8,0,0,0,0,1,1,1,0,2,11,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,0,0,7,11},
    (uint8_t []) {0,4,0,11,0,0,0,15,0,0,8,8,7,7,4,0,5,11,5,0,0,1,0,0,6,1,6,0,7,1,7,0,8,1,8,0,1,1,1,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,11,0,12,1,12,0,13,1,13,0,3,1,3,0,14,1,4,0,14,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,0,1,1,0,0,2,48,5,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,1,5,0,1,1,2,0,1,2,48,5,3},
    (uint8_t []) {0,2,0,10,0,0,0,12,0,0,7,7,2,0,3,11,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,6,0,1,1,1,0,7,1,7,0,8,1,8,0,9,1,9,0,10,11,10,0,11,11,2,0,11,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,16,16,2,0,3,1,3,0,0,1,0,0,1,1,1,0,4,1,2,0,4,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,11,2,0,3,11,3,0,4,1,4,0,0,1,0,0,5,1,1,0,5,1},
    (uint8_t []) {0,2,0,5,0,0,0,7,0,0,7,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,1,0,6,1,2,0,6,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,1,3,0,0,11,0,0,4,11,4,0,1,11,1,0,5,11,2,0,5,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,7,7,2,0,3,2,3,0,0,1,0,0,1,1,1,0,4,1,2,0,4,1,48,0,3},
    (uint8_t []) {0,4,0,0,0,0,0,4,0,0,16,26,16,26,0,0,1,1,1,0,2,1,2,0,3,1,0,0,3,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,8,8,2,0,3,1,3,0,4,1,4,0,0,1,0,0,1,1,1,0,5,1,2,0,5,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,8,8,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,0,0,7,0,0,0,7,0,0,0,0,1,11,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,11,5,0,6,1,0,0,6,1},
    (uint8_t []) {0,0,0,7,0,0,0,7,0,0,0,0,1,11,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,0,0,6,1},
    (uint8_t []) {0,6,0,13,0,0,0,19,0,0,8,7,7,7,7,7,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,2,1,2,0,10,1,10,0,11,1,11,0,3,1,3,0,12,1,12,0,13,1,13,0,4,1,4,0,14,1,14,0,15,1,15,0,5,1,5,0,16,1,16,0,17,1,17,0,0,1,0,0,18,1,6,0,18,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,1,0,5,11,2,0,5,11,0,1,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,0,1,1,0,0,11},
    (uint8_t []) {0,1,0,15,0,0,0,16,0,3,8,1,0,2,1,2,0,3,2,3,0,4,1,4,0,5,1,5,0,0,1,0,0,6,1,6,0,7,1,7,0,8,1,8,0,9,1,9,0,10,1,10,0,11,1,11,0,12,1,12,0,13,1,13,0,14,1,14,0,15,1,1,0,15,2,48,0,3,48,1,3,48,15,3},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,16,7,7,3,0,4,11,4,0,1,11,1,0,0,11,0,0,2,11,3,0,2,11},
    (uint8_t []) {0,1,0,15,0,0,0,16,0,4,8,1,0,2,1,2,0,3,2,3,0,4,1,4,0,5,1,5,0,6,1,6,0,0,1,0,0,7,1,7,0,8,1,8,0,9,1,9,0,10,1,10,0,11,1,11,0,12,2,12,0,13,1,13,0,14,1,14,0,15,1,1,0,15,2,48,0,3,48,1,3,48,11,3,48,15,3},
    (uint8_t []) {0,1,0,6,0,0,0,7,0,1,8,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,6,1,1,0,6,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,3,7,1,0,2,2,2,0,0,1,0,0,3,2,3,0,4,1,4,0,5,1,1,0,5,1,48,0,3,48,1,3,48,2,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,8,7,2,0,3,1,3,0,4,1,4,0,0,1,0,0,1,1,2,0,1,2,48,4,3},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,0,0,0,1,11,1,0,2,11,2,0,3,11,3,0,4,11,0,0,4,11},
    (uint8_t []) {0,0,0,8,0,0,0,8,0,0,0,0,1,1,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,0,0,7,1},
    (uint8_t []) {0,3,0,9,0,0,0,12,0,0,8,7,7,3,0,4,11,4,0,1,11,1,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,8,0,2,1,2,0,9,1,9,0,10,1,10,0,0,1,0,0,11,1,3,0,11,1},
    (uint8_t []) {0,1,0,3,0,0,0,4,0,0,8,1,0,2,1,2,0,0,1,0,0,3,1,1,0,3,1},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,3,26,7,7,3,0,4,2,4,0,5,1,5,0,1,2,1,0,0,1,0,0,2,1,3,0,2,11,0,0,254,0,1,1,48,2,3},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,5,27,7,7,3,0,4,2,4,0,1,1,1,0,0,1,0,0,2,1,2,0,5,2,3,0,5,1,0,0,254,0,2,1,48,0,3,48,4,3,48,5,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,1,3,0,0,11,0,0,4,11,4,0,1,11,1,0,5,11,2,0,5,11},
    (uint8_t []) {0,0,0,10,0,0,0,10,0,2,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,2,7,0,8,1,8,0,9,1,0,0,9,1,48,0,3,48,6,3},
    (uint8_t []) {0,1,0,13,0,0,0,14,0,0,8,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,0,1,0,0,8,1,8,0,9,1,9,0,10,1,10,0,11,1,11,0,12,1,12,0,13,1,1,0,13,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,11,5,0,1,1,2,0,1,2,48,5,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,3,2,3,0,0,1,0,0,4,1,4,0,1,1,1,0,5,1,2,0,5,1,48,0,3},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,2,7,1,0,2,1,2,0,3,1,3,0,4,11,4,0,0,1,1,0,0,2,0,0,1,48,4,3},
    (uint8_t []) {0,1,0,2,0,0,0,3,0,0,7,1,0,2,1,2,0,0,1,1,0,0,1},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,3,26,7,7,3,0,4,2,4,0,1,1,1,0,0,1,0,0,2,1,2,0,5,11,3,0,5,1,0,0,254,0,1,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,8,1,0,2,2,2,0,0,1,0,0,3,1,3,0,4,11,4,0,5,1,1,0,5,1,48,0,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,7,2,0,0,2,0,0,3,1,3,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1,48,0,3},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,1,7,1,0,2,1,2,0,3,1,3,0,4,11,4,0,0,1,1,0,0,2,48,4,3},
    (uint8_t []) {0,3,0,11,0,0,0,14,0,1,8,7,7,3,0,4,2,4,0,1,1,1,0,5,1,5,0,6,1,6,0,2,1,2,0,7,1,7,0,8,1,8,0,9,1,9,0,0,1,0,0,10,1,10,0,11,11,11,0,12,11,12,0,13,11,3,0,13,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,11,1,0,5,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,2,8,1,0,2,2,2,0,0,1,0,0,3,1,3,0,4,2,4,0,5,1,1,0,5,1,48,0,3,48,3,3},
    (uint8_t []) {0,4,0,12,0,0,0,16,0,0,7,7,7,7,4,0,5,11,5,0,0,11,0,0,6,11,6,0,7,11,7,0,8,11,8,0,1,11,1,0,9,11,9,0,10,11,10,0,11,11,11,0,2,11,2,0,12,11,12,0,13,11,13,0,14,11,14,0,3,11,3,0,15,11,4,0,15,11},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,1,1,1,0,5,1,2,0,5,1},
    (uint8_t []) {0,0,0,10,0,0,0,10,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,8,0,9,1,0,0,9,1,48,0,3},
    (uint8_t []) {0,1,0,6,0,0,0,7,0,0,7,1,0,2,11,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,1,0,6,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,16,2,0,3,11,3,0,0,11,0,0,1,11,1,0,4,11,2,0,4,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,8,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,11,4,0,5,1,1,0,5,1},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,0,0,0,1,1,1,0,2,11,2,0,3,1,3,0,4,1,0,0,4,11},
    (uint8_t []) {0,4,0,4,0,0,0,8,0,0,7,7,7,7,4,0,0,1,0,0,5,1,5,0,1,1,1,0,6,1,6,0,2,1,2,0,7,1,7,0,3,1,4,0,3,1},
    (uint8_t []) {0,3,0,13,0,0,0,16,0,0,8,7,7,3,0,4,11,4,0,5,11,5,0,0,1,0,0,6,1,6,0,7,11,7,0,8,11,8,0,9,11,9,0,10,1,10,0,11,1,11,0,12,1,12,0,1,1,1,0,13,1,13,0,14,1,14,0,2,1,2,0,15,1,3,0,15,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,16,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11,0,1,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,8,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,15,8,8,3,0,4,1,4,0,1,1,1,0,0,1,0,0,2,1,2,0,5,1,3,0,5,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,16,1,0,2,1,2,0,3,1,3,0,0,1,0,0,4,1,1,0,4,1},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,7,7,7,3,0,4,11,4,0,0,11,0,0,1,11,1,0,5,11,5,0,2,11,3,0,2,11},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,8,7,7,3,0,4,11,4,0,1,11,1,0,0,11,0,0,2,11,3,0,2,11},
    (uint8_t []) {0,3,0,13,0,0,0,16,0,0,8,7,7,3,0,4,11,4,0,5,11,5,0,0,1,0,0,6,1,6,0,7,11,7,0,8,11,8,0,9,11,9,0,10,1,10,0,11,1,11,0,1,1,1,0,12,1,12,0,13,1,13,0,2,1,2,0,14,1,14,0,15,1,3,0,15,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,1,8,1,0,2,2,2,0,0,1,0,0,3,1,3,0,4,1,1,0,4,1,48,0,3},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,1,7,1,0,2,2,2,0,0,1,0,0,3,1,3,0,4,1,1,0,4,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,8,1,0,2,1,2,0,0,11,0,0,3,11,3,0,4,11,4,0,5,11,1,0,5,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,2,2,0,0,1,0,0,3,1,3,0,4,1,4,0,5,1,1,0,5,1,48,0,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,7,8,2,0,0,2,0,0,3,1,3,0,4,1,4,0,1,1,2,0,1,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,11,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,1,0,0,1},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,0,0,0,1,11,1,0,2,11,2,0,3,1,3,0,4,1,4,0,5,1,0,0,5,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11,0,0,1},
    (uint8_t []) {0,1,0,2,0,0,0,3,0,0,26,1,0,2,1,2,0,0,1,1,0,0,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,16,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,11,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,8,7,7,3,0,1,11,1,0,4,11,4,0,0,11,0,0,2,11,3,0,2,11},
    (uint8_t []) {0,4,0,12,0,0,0,16,0,0,8,8,7,7,4,0,5,11,5,0,0,1,0,0,6,1,6,0,7,1,7,0,8,1,8,0,1,1,1,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,11,0,12,1,12,0,13,1,13,0,14,1,14,0,3,1,3,0,15,1,4,0,15,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,7,7,2,0,3,1,3,0,4,1,4,0,0,1,0,0,1,1,2,0,1,2,48,4,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,1,0,5,11,2,0,5,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,2,7,1,0,2,2,2,0,0,1,0,0,3,1,3,0,4,2,4,0,5,1,1,0,5,1,48,0,3,48,3,3},
    (uint8_t []) {0,1,0,6,0,0,0,7,0,0,8,1,0,2,1,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,1,0,6,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,3,7,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,2,1,0,4,1,48,0,3,48,3,3,48,4,3},
    (uint8_t []) {0,4,0,12,0,0,0,16,0,0,8,8,7,7,4,0,5,11,5,0,0,1,0,0,6,1,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,11,0,12,1,12,0,13,1,13,0,14,1,14,0,3,1,3,0,15,1,4,0,15,1},
    (uint8_t []) {0,0,0,7,0,0,0,7,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,0,0,6,1,48,0,3},
    (uint8_t []) {0,1,0,6,0,0,0,7,0,0,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,1,0,6,1},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,0,0,0,1,1,1,0,2,11,2,0,3,1,3,0,4,1,4,0,5,1,0,0,5,11},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,0,0,0,1,1,1,0,2,11,2,0,3,11,3,0,4,11,4,0,5,11,0,0,5,11},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,3,12,7,7,3,0,4,2,4,0,5,1,5,0,1,2,1,0,0,1,0,0,2,1,3,0,2,11,0,0,254,0,1,1,48,2,3},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,3,12,7,7,3,0,4,2,4,0,1,1,1,0,0,1,0,0,2,1,2,0,5,11,3,0,5,1,0,0,254,0,1,1,48,0,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,8,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,4,7,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,2,1,0,4,1,0,0,1,48,0,3,48,3,3,48,4,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,7,7,2,0,0,2,0,0,3,1,3,0,4,1,4,0,1,1,2,0,1,1,48,0,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,8,8,2,0,3,1,3,0,0,1,0,0,4,1,4,0,1,1,1,0,5,1,2,0,5,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,1,7,16,2,0,0,2,0,0,3,1,3,0,4,1,4,0,1,1,2,0,1,1,48,0,3},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,8,1,0,2,11,2,0,3,1,3,0,0,1,0,0,4,1,1,0,4,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,1,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,0,1,1,0,0,2,48,4,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,1,7,16,2,0,3,2,3,0,0,1,0,0,4,1,4,0,1,1,1,0,5,1,2,0,5,1,48,0,3},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,1,7,1,0,2,1,2,0,3,1,3,0,0,1,0,0,4,1,1,0,4,1,0,0,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,8,1,0,2,11,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,1,48,0,3},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,11,0,0,5,1,48,0,3},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,7,7,7,3,0,0,11,0,0,4,11,4,0,1,11,1,0,2,11,3,0,2,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,8,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,11,1,0,5,1,48,0,3},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,2,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,0,1,1,0,0,2,0,0,1,48,4,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,8,8,2,0,3,1,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,4,11,4,0,0,11,0,0,1,11,1,0,5,11,2,0,5,11},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,1,7,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,1,1,0,4,1,48,0,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,1,1,1,0,5,1,2,0,5,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,1,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,1,0,0,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,8,7,7,3,0,1,11,1,0,2,11,2,0,4,11,4,0,0,11,3,0,0,11},
    (uint8_t []) {0,3,0,3,0,0,0,6,0,0,7,7,7,3,0,0,11,0,0,4,11,4,0,1,11,1,0,5,11,5,0,2,11,3,0,2,11},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,11,4,0,5,1,0,0,5,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,8,1,0,2,11,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,1,0,5,11,0,0,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,11,3,0,0,11,0,0,1,11,1,0,4,11,2,0,4,11},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,16,7,7,3,0,1,11,1,0,2,11,2,0,4,11,4,0,0,11,3,0,0,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,8,1,0,2,2,2,0,0,1,0,0,3,1,3,0,4,1,4,0,5,1,1,0,5,1,48,0,3},
    (uint8_t []) {0,0,0,7,0,0,0,7,0,0,0,0,1,1,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,5,0,6,1,0,0,6,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,8,8,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,11,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,1,0,5,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,8,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,1,2,0,3,1,3,0,4,1,4,0,0,1,1,0,0,11},
    (uint8_t []) {0,4,0,1,0,0,0,5,0,0,7,7,7,7,4,0,0,11,0,0,1,11,1,0,2,11,2,0,3,11,4,0,3,11},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,7,1,0,2,11,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,1,0,5,11,0,0,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,11,2,0,3,1,3,0,0,1,0,0,4,1,1,0,4,1},
    (uint8_t []) {0,2,0,6,0,0,0,8,0,0,7,7,2,0,3,1,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,6,0,1,1,1,0,7,1,2,0,7,1},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,0,0,0,1,11,1,0,2,1,2,0,3,1,3,0,4,1,0,0,4,1},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,0,0,0,1,11,1,0,2,1,2,0,3,1,3,0,4,11,4,0,5,1,0,0,5,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,1},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,3,0,0,1,1,1,0,2,2,2,0,3,1,3,0,4,1,4,0,5,1,0,0,5,2,48,0,3,48,1,3,48,5,3},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,16,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,1,8,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,1,1,0,4,1,48,0,3},
    (uint8_t []) {0,3,0,9,0,0,0,12,0,0,8,7,7,3,0,4,1,4,0,0,1,0,0,5,1,5,0,6,1,6,0,7,1,7,0,1,1,1,0,8,1,8,0,9,1,9,0,10,1,10,0,11,1,11,0,2,1,3,0,2,11},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,2,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,2,4,0,5,1,0,0,5,1,48,0,3,48,3,3},
    (uint8_t []) {0,3,0,10,0,0,0,13,0,0,8,7,7,3,0,4,11,4,0,1,11,1,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,8,0,2,1,2,0,9,1,9,0,10,1,10,0,11,1,11,0,0,1,0,0,12,1,3,0,12,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,7,7,7,3,0,0,11,0,0,4,11,4,0,1,11,1,0,2,11,3,0,2,11},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,1,1,0,4,1},
    (uint8_t []) {0,0,0,4,0,0,0,4,0,0,0,0,1,1,1,0,2,1,2,0,3,1,0,0,3,1},
    (uint8_t []) {0,2,0,6,0,0,0,8,0,0,7,8,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,6,1,6,0,1,1,1,0,7,1,2,0,7,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,5,0,1,11,2,0,1,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,1,4,0,5,1,1,0,5,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,0,0,4,1,48,0,3},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,8,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,1,4,0,5,1,1,0,5,1},
    (uint8_t []) {0,2,0,9,0,0,0,11,0,0,8,7,2,0,3,1,3,0,4,11,4,0,5,1,5,0,0,1,0,0,6,1,6,0,7,1,7,0,8,1,8,0,1,1,1,0,9,1,9,0,10,1,2,0,10,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,11,3,0,0,11,0,0,1,11,1,0,4,11,2,0,4,11},
    (uint8_t []) {0,1,0,2,0,0,0,3,0,0,8,1,0,2,1,2,0,0,1,1,0,0,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,8,2,0,3,11,3,0,0,11,0,0,1,11,1,0,4,11,2,0,4,11},
    (uint8_t []) {0,2,0,7,0,0,0,9,0,0,8,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,5,0,6,1,6,0,1,1,1,0,7,1,7,0,8,1,2,0,8,1},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,0,0,0,1,11,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,0,0,5,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,1,8,1,0,2,2,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,1,48,0,3},
    (uint8_t []) {0,3,0,5,0,0,0,8,0,0,16,7,8,3,0,4,11,4,0,0,1,0,0,1,1,1,0,5,1,5,0,6,1,6,0,7,1,7,0,2,1,3,0,2,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,8,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,16,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,11,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,1},
    (uint8_t []) {0,3,0,2,0,0,0,5,0,0,7,7,7,3,0,4,11,4,0,0,11,0,0,1,11,1,0,2,11,3,0,2,11},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,8,1,0,2,11,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,1,0,5,11},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,16,1,0,2,11,2,0,3,11,3,0,0,11,0,0,4,11,1,0,4,11},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,8,1,0,2,11,2,0,3,11,3,0,0,11,0,0,4,11,1,0,4,11},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,1,0,5,11,2,0,5,11},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,11,2,0,3,11,3,0,0,11,0,0,4,11,1,0,4,11},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,8,8,2,0,3,11,3,0,0,1,0,0,4,1,4,0,1,1,2,0,1,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,11,2,0,3,11,3,0,0,11,0,0,4,11,1,0,4,11},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,8,1,0,2,11,2,0,0,1,0,0,3,1,3,0,4,1,1,0,4,1},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,5,0,1,1,2,0,1,1},
    (uint8_t []) {0,3,0,11,0,0,0,14,0,0,7,8,8,3,0,4,11,4,0,1,1,1,0,5,1,5,0,6,1,6,0,7,1,7,0,8,1,8,0,9,1,9,0,2,1,2,0,10,1,10,0,11,1,11,0,12,1,12,0,0,1,0,0,13,1,3,0,13,1},
    (uint8_t []) {0,0,0,3,0,0,0,3,0,0,0,0,1,1,1,0,2,1,0,0,2,1},
    (uint8_t []) {0,3,0,9,0,0,0,12,0,0,8,7,7,3,0,4,11,4,0,0,1,0,0,5,1,5,0,6,1,6,0,1,1,1,0,7,1,7,0,8,1,8,0,9,1,9,0,10,1,10,0,2,1,2,0,11,1,3,0,11,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,7,1,0,2,1,2,0,3,1,3,0,0,1,0,0,4,1,1,0,4,1},
    (uint8_t []) {0,1,0,3,0,0,0,4,0,0,7,1,0,2,1,2,0,0,1,0,0,3,1,1,0,3,1},
    (uint8_t []) {0,2,0,6,0,0,0,8,0,0,8,7,2,0,3,11,3,0,0,1,0,0,4,1,4,0,5,1,5,0,6,1,6,0,1,1,1,0,7,1,2,0,7,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,1,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,1},
    (uint8_t []) {0,2,0,3,0,0,0,5,0,0,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,2,0,1,11},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,1,0,0,1,2,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,0,0,5,1,48,0,3},
    (uint8_t []) {0,2,0,4,0,0,0,6,0,0,7,7,2,0,3,11,3,0,0,11,0,0,4,11,4,0,1,11,1,0,5,11,2,0,5,11},
    (uint8_t []) {0,0,0,5,0,0,0,5,0,0,0,0,1,1,1,0,2,1,2,0,3,1,3,0,4,1,0,0,4,1},
    (uint8_t []) {0,1,0,4,0,0,0,5,0,0,8,1,0,2,1,2,0,3,1,3,0,0,1,0,0,4,1,1,0,4,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,7,1,0,2,11,2,0,3,11,3,0,0,11,0,0,4,11,4,0,5,11,1,0,5,11},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,0,0,0,1,1,1,0,2,1,2,0,3,1,3,0,4,1,4,0,5,1,0,0,5,1},
    (uint8_t []) {0,1,0,5,0,0,0,6,0,0,8,1,0,2,1,2,0,3,1,3,0,0,1,0,0,4,1,4,0,5,1,1,0,5,1},
    (uint8_t []) {0,0,0,6,0,0,0,6,0,0,0,0,1,11,1,0,2,11,2,0,3,11,3,0,4,11,4,0,5,11,0,0,5,11}
};


void crng_fingerprint_init(void)
{
    PG_MEMCONTEXT_BEGIN(TopMemoryContext);
    for(int i = 0; i < PATTERN_COUNT; i++)
        molecule_simple_init(patternMolecule + i, patterns[i]);
    PG_MEMCONTEXT_END();
}


std::map<uint32_t, int> crng_fingerprint_get(const Molecule *molecule, BitInfo *info)
{
    if(unlikely(initialized == false))
    {
        crng_fingerprint_init();
        initialized = true;
    }


    std::map<uint32_t, int> fp;

    for(int i = 0; i < PATTERN_COUNT; i++)
    {
        SubstructureMatch substructure(patternMolecule + i);
        std::vector<std::vector<int>> matches = substructure.match(molecule, 256);

        if(matches.size())
        {
            fp[i] = matches.size();

            if(info)
            {
                for(auto &match : matches)
                    (*info)[i].insert(match.begin(), match.end());
            }
        }
    }

    return fp;
}
