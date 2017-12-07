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
#endif

#define PG_MEMCONTEXT_BEGIN(context)    do { MemoryContext old = MemoryContextSwitchTo(context)
#define PG_MEMCONTEXT_END()             MemoryContextSwitchTo(old);} while(0)


inline void createBasePath(void)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));

    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);

    char *path = (char *) palloc(basePathLength +  databaseLength + 2);
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


inline char *getFilePath(char* fileName)
{
    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));

    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);
    size_t fileNameLength = strlen(fileName);

    char *path = (char *) palloc(basePathLength +  databaseLength + fileNameLength + 3);
    char *data = path;

    memcpy(data, DataDir, basePathLength);
    data += basePathLength;

    *(data++) = '/';

    memcpy(data, database->data, databaseLength);
    data += databaseLength;

    *(data++) = '/';

    memcpy(data, fileName, fileNameLength);
    data += fileNameLength;

    *data = '\0';

    return path;
}

#endif /* SACHEM_H_ */
