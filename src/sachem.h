#ifndef SACHEM_H_
#define SACHEM_H_

#include <postgres.h>
#include <fmgr.h>
#include <utils/builtins.h>
#include <miscadmin.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>


#define INDEX_PREFIX_NAME       "lucene"
#define COMPOUNDS_TABLE         "sachem.compounds"
#define MOLECULE_ERRORS_TABLE   "sachem.sachem_molecule_errors"
#define AUDIT_TABLE             "sachem.sachem_compound_audit"
#define INDEX_TABLE             "sachem.sachem_index"


#if PG_VERSION_NUM < 100000
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)
#define shm_toc_lookup_key(toc,key) shm_toc_lookup((toc),(key))
#elif PG_VERSION_NUM < 110000
#define CreateParallelContextForExternalFunction(contex,library,function) CreateParallelContext((contex),(library),(function))
#define shm_toc_lookup_key(toc,key) shm_toc_lookup((toc),(key),false)
#else
#define CreateParallelContextForExternalFunction(contex,library,function) CreateParallelContext((contex),(library),(function),false)
#define shm_toc_lookup_key(toc,key) shm_toc_lookup((toc),(key),false)
#endif

#define PG_MEMCONTEXT_BEGIN(context)    do { MemoryContext old = MemoryContextSwitchTo(context)
#define PG_MEMCONTEXT_END()             MemoryContextSwitchTo(old);} while(0)


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


static inline char *get_index_name(const char *prefix, int indexNumber)
{
    char *indexPath = (char *) palloc(strlen(prefix) + 24);
    sprintf(indexPath, "%s-%i", prefix, indexNumber);

    return indexPath;
}


static inline char *get_index_path(const char *prefix, int indexNumber)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));
    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);
    size_t prefixLength = strlen(prefix);

    char *indexPath = (char *) palloc(basePathLength + databaseLength + prefixLength + 24);
    sprintf(indexPath, "%s/%s/%s-%i", DataDir, database->data, prefix, indexNumber);

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
