#ifndef STATS_H__
#define STATS_H__

#include "molecule.h"


typedef void *Stats;


Stats *stats_create();
void stats_delete(Stats *stats);
bool stats_add(Stats *stats, const Molecule *molecule);
bool stats_merge(Stats *stats, Stats *substats);
bool stats_write(Stats *stats, const char *name, int limit);

#endif /*STATS_H__*/
