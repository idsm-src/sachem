#ifndef SACHEM_H_
#define SACHEM_H_

#include <postgres.h>
#include <fmgr.h>
#include <utils/builtins.h>
#include <miscadmin.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>


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


static inline Datum SPI_get_value(HeapTuple row, TupleDesc rowdesc, int colnumber)
{
    bool isNullFlag;

    Datum result = SPI_getbinval(row, rowdesc, colnumber, &isNullFlag);

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);

    return result;
}


static inline void create_base_directory(char *indexName)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));

    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);
    size_t indexNameLength = strlen(indexName);

    char *fragment[] = { "sachem", database->data, indexName };
    char *path = (char *) palloc(basePathLength + databaseLength + indexNameLength + 10);
    char *data = path;

    memcpy(data, DataDir, basePathLength);
    data += basePathLength;

    for(int i = 0; i < 3; i++)
    {
        *(data++) = '/';

        size_t length = strlen(fragment[i]);

        memcpy(data, fragment[i], length + 1);
        data += length;

        if(unlikely(mkdir(path, S_IRUSR | S_IWUSR | S_IXUSR) == -1 && errno != EEXIST))
            elog(ERROR, "%s: mkdir() failed", __func__);
    }
}


static inline bool is_index_name(const char *name)
{
    return !strncmp(name, "lucene-", 7);
}


static inline char *get_index_name(int indexNumber)
{
    char *indexPath = (char *) palloc(30);
    sprintf(indexPath, "lucene-%i", indexNumber);

    return indexPath;
}


static inline char *get_index_path(const char *indexName, int indexNumber)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));
    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);
    size_t indexNameLength = strlen(indexName);

    char *indexPath = (char *) palloc(basePathLength + databaseLength + indexNameLength + 38);
    sprintf(indexPath, "%s/sachem/%s/%s/lucene-%i", DataDir, database->data, indexName, indexNumber);

    return indexPath;
}


static inline char *get_base_path(const char *indexName)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));
    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);
    size_t indexNameLength = strlen(indexName);

    char *path = (char *) palloc(basePathLength + databaseLength + indexNameLength + 10);
    sprintf(path, "%s/sachem/%s/%s", DataDir, database->data, indexName);

    return path;
}

#endif /* SACHEM_H_ */
