#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <utils/builtins.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include "java.h"
#include "sachem.h"


#define SYNC_FETCH_SIZE         100000


static bool initialized = false;
static jclass indexerClass;
static jmethodID constructor;
static jmethodID beginMethod;
static jmethodID addMethod;
static jmethodID deleteMethod;
static jmethodID optimizeMethod;
static jmethodID commitMethod;
static jmethodID rollbackMethod;


static void indexer_java_init()
{
    if(likely(initialized))
        return;

    java_init();


    indexerClass = (*env)->FindClass(env, "cz/iocb/sachem/lucene/Indexer");
    java_check_exception(__func__);

    constructor = (*env)->GetMethodID(env, indexerClass, "<init>", "()V");
    java_check_exception(__func__);

    beginMethod = (*env)->GetMethodID(env, indexerClass, "begin", "(Ljava/lang/String;IID)V");
    java_check_exception(__func__);

    addMethod = (*env)->GetMethodID(env, indexerClass, "add", "(I[B)Ljava/lang/String;");
    java_check_exception(__func__);

    deleteMethod = (*env)->GetMethodID(env, indexerClass, "delete", "(I)V");
    java_check_exception(__func__);

    optimizeMethod = (*env)->GetMethodID(env, indexerClass, "optimize", "()V");
    java_check_exception(__func__);

    commitMethod = (*env)->GetMethodID(env, indexerClass, "commit", "()V");
    java_check_exception(__func__);

    rollbackMethod = (*env)->GetMethodID(env, indexerClass, "rollback", "()V");
    java_check_exception(__func__);

    initialized = true;
}


static jobject indexer_init()
{
    indexer_java_init();

    jobject indexer = (*env)->NewObject(env, indexerClass, constructor);
    java_check_exception(__func__);

    return indexer;
}


static void indexer_terminate(jobject indexer)
{
    JavaDeleteRef(indexer);
}


static void indexer_begin(jobject indexer, const char *path, int segments, int bufferedDocs, double bufferSize)
{
    jstring folder = NULL;

    PG_TRY();
    {
        folder = (*env)->NewStringUTF(env, path);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, indexer, beginMethod, folder, segments, bufferedDocs, bufferSize);
        java_check_exception(__func__);

        JavaDeleteRef(folder);
    }
    PG_CATCH();
    {
        JavaDeleteRef(folder);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


static VarChar *indexer_add(jobject indexer, int32_t id, VarChar *molfile)
{
    jbyteArray molfileArray = NULL;
    jstring message = NULL;
    VarChar *result = NULL;

    PG_TRY();
    {
        size_t length = VARSIZE(molfile) - VARHDRSZ;

        molfileArray = (jbyteArray) (*env)->NewByteArray(env, length);
        java_check_exception(__func__);

        (*env)->SetByteArrayRegion(env, molfileArray, 0, length, (jbyte*) VARDATA(molfile));
        java_check_exception(__func__);

        message = (jstring)(*env)->CallObjectMethod(env, indexer, addMethod, (jint) id, molfileArray);
        java_check_exception(__func__);

        if(unlikely(message != NULL))
        {
            jsize length = (*env)->GetStringUTFLength(env, message);

            result = (VarChar *) palloc(VARHDRSZ + length);
            SET_VARSIZE(result, VARHDRSZ + length);

            const char *mstr = (*env)->GetStringUTFChars(env, message, NULL);

            if(unlikely(mstr == NULL))
                elog(ERROR, "unknown jvm error");

            memcpy(VARDATA(result), mstr, length);

            (*env)->ReleaseStringUTFChars(env, message, mstr);
            JavaDeleteRef(message);
        }

        JavaDeleteRef(molfileArray);
    }
    PG_CATCH();
    {
        JavaDeleteRef(molfileArray);
        JavaDeleteRef(message);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return result;
}


static void indexer_delete(jobject indexer, int32_t id)
{
    (*env)->CallVoidMethod(env, indexer, deleteMethod, (jint) id);
    java_check_exception(__func__);
}


static void indexer_optimize(jobject indexer)
{
    (*env)->CallVoidMethod(env, indexer, optimizeMethod);
    java_check_exception(__func__);
}


static void indexer_commit(jobject indexer)
{
    (*env)->CallVoidMethod(env, indexer, commitMethod);
    java_check_exception(__func__);
}


static void indexer_rollback(jobject indexer)
{
    (*env)->CallVoidMethod(env, indexer, rollbackMethod);
    java_check_exception(__func__);
}


static void link_directory(const char *oldPath, const char *newPath)
{
    int newdirfd = -1;
    int olddirfd = -1;
    DIR *dp = NULL;


    PG_TRY();
    {
        if(unlikely(mkdir(newPath, 0700) != 0))
            elog(ERROR, "%s: mkdir() failed", __func__);


        if(likely(oldPath != NULL))
        {
            if(unlikely((newdirfd = open(newPath, O_DIRECTORY)) == -1))
                elog(ERROR, "%s: open() failed", __func__);

            if(unlikely((olddirfd = open(oldPath, O_DIRECTORY)) == -1))
                elog(ERROR, "%s: open() failed", __func__);

            if(unlikely((dp = fdopendir(olddirfd)) == NULL))
                elog(ERROR, "%s: fdopendir() failed", __func__);


            struct dirent *ep;

            while((ep = readdir(dp)))
            {
                if(unlikely(ep->d_name[0] == '.'))
                    continue;

                if(unlikely(linkat(olddirfd, ep->d_name, newdirfd, ep->d_name, 0) != 0))
                    elog(ERROR, "%s: linkat() failed", __func__);
            }

            closedir(dp);
            dp = NULL;

            close(olddirfd);
            olddirfd = -1;

            close(newdirfd);
            newdirfd = -1;
        }
    }
    PG_CATCH();
    {
        if(dp != NULL)
            closedir(dp);

        if(olddirfd != -1)
            close(olddirfd);

        if(newdirfd != -1)
            close(newdirfd);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


static void delete_directory(int parent, const char *path)
{
    int dirfd = -1;
    DIR *dp = NULL;

    PG_TRY();
    {
        if(unlikely((dirfd = openat(parent, path, O_DIRECTORY)) == -1))
            elog(ERROR, "%s: openat() failed", __func__);

        if(unlikely((dp = fdopendir(dirfd)) == NULL))
            elog(ERROR, "%s: fdopendir() failed", __func__);


        struct dirent *ep;

        while((ep = readdir(dp)))
        {
            if(unlikely(ep->d_name[0] == '.'))
                continue;

            if(unlikely(unlinkat(dirfd, ep->d_name, 0) != 0))
                elog(ERROR, "%s: unlinkat() failed", __func__);
        }

        closedir(dp);
        dp = NULL;

        close(dirfd);
        dirfd = -1;

        if(unlikely(unlinkat(parent, path, AT_REMOVEDIR) != 0))
           elog(ERROR, "%s: unlinkat() failed", __func__);
    }
    PG_CATCH();
    {
        if(dp != NULL)
            closedir(dp);

        if(dirfd != -1)
            close(dirfd);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


PG_FUNCTION_INFO_V1(sync_data);
Datum sync_data(PG_FUNCTION_ARGS)
{
    VarChar *index = PG_GETARG_VARCHAR_P(0);
    bool verbose = PG_GETARG_BOOL(1);
    bool optimize = PG_GETARG_BOOL(2);


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);


    /* load configuration */
    if(unlikely(SPI_execute_with_args("select id, version, quote_ident(schema_name), quote_ident(table_name), "
            "quote_ident(id_column), quote_ident(molfile_column), segments, buffered_docs, buffer_size from "
            "sachem.configuration where index_name = $1", 1,
            (Oid[]) { VARCHAROID }, (Datum[]) { PointerGetDatum(index) }, NULL, true, 1) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 9))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    Datum indexId = SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1);
    Datum version = SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2);
    char *schemaName = text_to_cstring(DatumGetVarCharP(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3)));
    char *tableName = text_to_cstring(DatumGetVarCharP(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 4)));
    char *idColumn = text_to_cstring(DatumGetVarCharP(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 5)));
    char *molfileColumn = text_to_cstring(DatumGetVarCharP(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 6)));
    int32 segments = DatumGetInt32(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 7));
    int32 bufferedDocs = DatumGetInt32(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 8));
    float8 bufferSize = DatumGetFloat8(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 9));
    char *indexName = text_to_cstring(index);


    /* clone old index */
    int oldIndexVersion = DatumGetInt32(version);
    int indexVersion = oldIndexVersion + 1;

    char *indexPath = get_index_path(indexName, indexVersion);

    if(unlikely(oldIndexVersion == 0))
        create_base_directory(indexName);
    else
        link_directory(get_index_path(indexName, oldIndexVersion), indexPath);


    /* update version */
    if(unlikely(SPI_execute_with_args("update sachem.configuration set version = $2 where id = $1;", 2,
            (Oid[]) { INT4OID, INT4OID }, (Datum[]) { indexId, Int32GetDatum(indexVersion) }, NULL, false, 0) != SPI_OK_UPDATE))
        elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);



    jobject indexer = indexer_init();

    PG_TRY();
    {
        indexer_begin(indexer, indexPath, segments, bufferedDocs, bufferSize);

        /* delete unnecessary data */
        Portal auditCursor = SPI_cursor_open_with_args(NULL,
                "select id from sachem.compound_audit where not stored and index = $1",
                1, (Oid[]) { INT4OID }, (Datum[]) { indexId }, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

        while(true)
        {
            SPI_cursor_fetch(auditCursor, true, SYNC_FETCH_SIZE);

            if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
                elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

            if(unlikely(SPI_processed == 0))
                break;

            for(size_t i = 0; i < SPI_processed; i++)
            {
                CHECK_FOR_INTERRUPTS();

                indexer_delete(indexer, DatumGetInt32(SPI_get_value(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1)));
            }
        }

        SPI_cursor_close(auditCursor);


        /* convert new data */
        char *query = (char *) palloc(110 + 2 * strlen(idColumn) + strlen(molfileColumn) + strlen(schemaName) + strlen(tableName));
        sprintf(query, "select cmp.%s, cmp.%s from %s.%s cmp, sachem.compound_audit aud where cmp.%s = aud.id and "
                "aud.stored and aud.index = $1", idColumn, molfileColumn, schemaName, tableName, idColumn);

        Portal compoundCursor = SPI_cursor_open_with_args(NULL, query, 1, (Oid[]) { INT4OID }, (Datum[]) { indexId },
                NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

        uint64 count = 0;

        while(true)
        {
            SPI_cursor_fetch(compoundCursor, true, SYNC_FETCH_SIZE);

            if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
                elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

            if(unlikely(SPI_processed == 0))
                break;

            uint64 processed = SPI_processed;
            SPITupleTable *tuptable = SPI_tuptable;

            for(int i = 0; i < processed; i++)
            {
                CHECK_FOR_INTERRUPTS();

                Datum id = SPI_get_value(tuptable->vals[i], tuptable->tupdesc, 1);
                Datum molfile = SPI_get_value(tuptable->vals[i], tuptable->tupdesc, 2);

                VarChar *message = indexer_add(indexer, DatumGetInt32(id), DatumGetVarCharP(molfile));

                if(unlikely(message != NULL))
                {
                    char *notice = text_to_cstring(message);
                    elog(NOTICE, "%i: %s", DatumGetInt32(id), notice);
                    pfree(notice);

                    if(unlikely(SPI_execute_with_args(
                            "insert into sachem.compound_errors (index, compound, message) values ($1,$2,$3)", 3,
                            (Oid[]) { INT4OID, INT4OID, TEXTOID }, (Datum[]) { indexId, id, PointerGetDatum(message) },
                            NULL, false, 0) != SPI_OK_INSERT))
                        elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);

                    pfree(message);
                }
            }

            SPI_freetuptable(tuptable);

            count += processed;

            if(verbose)
                elog(NOTICE, "already processed: %lu", count);
        }

        SPI_cursor_close(compoundCursor);

        if(unlikely(SPI_execute_with_args("delete from sachem.compound_audit where index = $1", 1, (Oid[]) { INT4OID },
                (Datum[]) { indexId }, NULL, false, 0) != SPI_OK_DELETE))
            elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);



        if(optimize)
            indexer_optimize(indexer);

        indexer_commit(indexer);
    }
    PG_CATCH();
    {
        indexer_rollback(indexer);
        indexer_terminate(indexer);

        PG_RE_THROW();
    }
    PG_END_TRY();

    indexer_terminate(indexer);

    SPI_finish();
    PG_RETURN_VOID();
}


PG_FUNCTION_INFO_V1(cleanup);
Datum cleanup(PG_FUNCTION_ARGS)
{
    VarChar *index = PG_GETARG_VARCHAR_P(0);
    char *indexName = text_to_cstring(index);


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    int version = 0;

    if(unlikely(SPI_execute_with_args("select version from sachem.configuration where index_name = $1",
            1, (Oid[]) { VARCHAROID }, (Datum[]) { PointerGetDatum(index) }, NULL, false, 1) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);

    if(likely(SPI_processed != 0))
    {
        if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
            elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

        version = DatumGetInt32(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1));
    }

    SPI_finish();


    char *luceneIndexName = get_index_name(version);
    int dirfd = -1;
    DIR *dp = NULL;

    if(unlikely((dirfd = open(get_base_path(indexName), O_DIRECTORY)) == -1 && errno != ENOENT))
        elog(ERROR, "%s: open() failed", __func__);

    if(dirfd == -1)
        PG_RETURN_VOID();


    PG_TRY();
    {
        if(unlikely((dp = fdopendir(dirfd)) == NULL))
            elog(ERROR, "%s: fdopendir() failed", __func__);


        struct dirent *ep;

        while((ep = readdir(dp)))
        {
            if(likely(is_index_name(ep->d_name)))
            {
                if(unlikely(strcmp(ep->d_name, luceneIndexName)))
                {
                    elog(NOTICE, "delete lucene index '%s'", ep->d_name);
                    delete_directory(dirfd, ep->d_name);
                }
            }
            else if(unlikely(strcmp(ep->d_name, ".") && strcmp(ep->d_name, "..")))
            {
                elog(WARNING, "unknown content '%s'", ep->d_name);
            }
        }

        closedir(dp);
        dp = NULL;

        close(dirfd);
        dirfd = -1;
    }
    PG_CATCH();
    {
        if(dp != NULL)
            closedir(dp);

        if(dirfd != -1)
            close(dirfd);

        PG_RE_THROW();
    }
    PG_END_TRY();

    PG_RETURN_VOID();
}
