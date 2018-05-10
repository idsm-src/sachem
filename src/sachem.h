#ifndef SACHEM_H_
#define SACHEM_H_

#include <postgres.h>
#include <fmgr.h>
#include <utils/builtins.h>
#include <miscadmin.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef __cplusplus
#define restrict
#endif

#if PG_VERSION_NUM < 100000
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)
#define shm_toc_lookup_key(toc,key) shm_toc_lookup((toc),(key))
#else
#define CreateParallelContextForExternalFunction(contex,library,function) CreateParallelContext((contex),(library),(function))
#define shm_toc_lookup_key(toc,key) shm_toc_lookup((toc),(key),false)
#endif

#define PG_MEMCONTEXT_BEGIN(context)    do { MemoryContext old = MemoryContextSwitchTo(context)
#define PG_MEMCONTEXT_END()             MemoryContextSwitchTo(old);} while(0)


#define SAFE_CPP_BEGIN                                                                              \
    bool badAlloc = false;                                                                          \
    try                                                                                             \
    {


#define SAFE_CPP_END                                                                                \
    }                                                                                               \
    catch(std::bad_alloc &e)                                                                        \
    {                                                                                               \
        badAlloc = true;                                                                            \
    }                                                                                               \
    catch(...)                                                                                      \
    {                                                                                               \
    }                                                                                               \
    if(badAlloc)                                                                                    \
        ereport(ERROR, (errcode(ERRCODE_OUT_OF_MEMORY), errmsg("out of memory")));                  \
    else                                                                                            \
        elog(ERROR, "%s: unexpected exception", __func__);


static inline void create_base_directory(void)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));

    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);

    char *path = (char *) palloc(basePathLength + databaseLength + 2);
    char *data = path;

    memcpy(data, DataDir, basePathLength);
    data += basePathLength;

    *(data++) = '/';

    memcpy(data, database->data, databaseLength);
    data += databaseLength;

    *data = '\0';

    if(mkdir(path, S_IRUSR | S_IWUSR | S_IXUSR) == -1 && errno != EEXIST)
        elog(ERROR, "%s: mkdir() failed", __func__);
}


static inline char *get_index_path(const char *prefix, const char *suffix, int indexNumber)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));
    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);
    size_t prefixLength = strlen(prefix);
    size_t suffixLength = strlen(suffix);

    char *indexPath = (char *) palloc(basePathLength + databaseLength + prefixLength + suffixLength + 24);
    sprintf(indexPath, "%s/%s/%s-%i%s", DataDir, database->data, prefix, indexNumber, suffix);

    return indexPath;
}


static inline char *get_subindex_path(const char *prefix, const char *suffix, int indexNumber, int subNumber)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));
    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);
    size_t prefixLength = strlen(prefix);
    size_t suffixLength = strlen(suffix);

    char *indexPath = (char *) palloc(basePathLength + databaseLength + prefixLength + suffixLength + 64);
    sprintf(indexPath, "%s/%s/%s-%i.%i%s", DataDir, database->data, prefix, indexNumber, subNumber, suffix);

    return indexPath;
}


static inline char *get_file_path(const char *name)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));
    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);
    size_t nameLength = strlen(name);

    char *filePath = (char *) palloc(basePathLength + databaseLength + nameLength + 3);
    sprintf(filePath, "%s/%s/%s", DataDir, database->data, name);

    return filePath;
}

#endif /* SACHEM_H_ */
