#ifndef PGCHEM_H_
#define PGCHEM_H_

#define SHOW_STATS              0
#define USE_COUNT_FINGERPRINT   1
#define BITSET_ASSERT           0

#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

#define PG_MEMCONTEXT_BEGIN(context)    do { MemoryContext old = MemoryContextSwitchTo(mcxt)
#define PG_MEMCONTEXT_END()             MemoryContextSwitchTo(old);} while(0)

#endif /* PGCHEM_H_ */
