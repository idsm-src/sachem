#ifndef STATS_H__
#define STATS_H__

#include "molecule.h"


typedef void *Stats;


typedef struct
{
    uint32_t fp;
    uint32_t count;
} StatItem;


Stats *stats_create();
void stats_delete(Stats *stats);
void stats_add(Stats *stats, const Molecule *molecule);
void stats_merge(Stats *stats, StatItem *items, size_t size);
size_t stats_get_items(Stats *stats, StatItem **items);
void stats_write(Stats *stats, const char *name, int limit);

#endif /*STATS_H__*/
