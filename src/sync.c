#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <unistd.h>
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

    beginMethod = (*env)->GetMethodID(env, indexerClass, "begin", "(Ljava/lang/String;)V");
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


static void indexer_begin(jobject indexer, const char *path)
{
    jstring folder = NULL;

    PG_TRY();
    {
        folder = (*env)->NewStringUTF(env, path);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, indexer, beginMethod, folder);
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

        if(message != NULL)
        {
            jsize length = (*env)->GetStringUTFLength(env, message);

            result = (VarChar *) palloc(VARHDRSZ + length);
            SET_VARSIZE(result, VARHDRSZ + length);

            const char *mstr = (*env)->GetStringUTFChars(env, message, NULL);

            if(mstr == NULL)
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
        if(mkdir(newPath, 0700) != 0)
            elog(ERROR, "%s: mkdir() failed", __func__);


        if(oldPath != NULL)
        {
            if((newdirfd = open(newPath, O_DIRECTORY)) == -1)
                elog(ERROR, "%s: open() failed", __func__);

            if((olddirfd = open(oldPath, O_DIRECTORY)) == -1)
                elog(ERROR, "%s: open() failed", __func__);

            if((dp = fdopendir(olddirfd)) == NULL)
                elog(ERROR, "%s: fdopendir() failed", __func__);


            struct dirent *ep;

            while((ep = readdir(dp)))
            {
                if(ep->d_name[0] == '.')
                    continue;

                if(linkat(olddirfd, ep->d_name, newdirfd, ep->d_name, 0) != 0)
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


static void delete_directory(const char *path)
{
    int dirfd = -1;
    DIR *dp = NULL;

    PG_TRY();
    {
        if((dirfd = open(path, O_DIRECTORY)) == -1)
            elog(ERROR, "%s: open() failed", __func__);

        if((dp = fdopendir(dirfd)) == NULL)
            elog(ERROR, "%s: fdopendir() failed", __func__);


        struct dirent *ep;

        while((ep = readdir(dp)))
        {
            if(ep->d_name[0] == '.')
                continue;

            if(unlinkat(dirfd, ep->d_name, 0) != 0)
                elog(ERROR, "%s: unlinkat() failed", __func__);
        }

        closedir(dp);
        dp = NULL;

        close(dirfd);
        dirfd = -1;

        if(rmdir(path) != 0)
           elog(ERROR, "%s: rmdir() failed", __func__);
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


PG_FUNCTION_INFO_V1(sachem_sync_data);
Datum sachem_sync_data(PG_FUNCTION_ARGS)
{
    create_base_directory();

    bool verbose = PG_GETARG_BOOL(0);
    bool optimize = PG_GETARG_BOOL(1);


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    bool isNullFlag;


    /* clone old index */
    int indexNumber = 0;
    char *oldIndexPath = NULL;

    if(unlikely(SPI_exec("select id from " INDEX_TABLE, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(SPI_processed != 0)
    {
        if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
            elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

        Datum number = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE) || isNullFlag)
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        oldIndexPath = get_index_path(INDEX_PREFIX_NAME, DatumGetInt32(number));
        indexNumber = DatumGetInt32(number) + 1;
    }

    char *indexPath = get_index_path(INDEX_PREFIX_NAME, indexNumber);

    link_directory(oldIndexPath, indexPath);

    if(unlikely(SPI_exec("delete from " INDEX_TABLE, 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(SPI_execute_with_args("insert into " INDEX_TABLE " (id) values ($1)", 1, (Oid[]) { INT4OID },
            (Datum[]) {Int32GetDatum(indexNumber)}, NULL, false, 0) != SPI_OK_INSERT)
        elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);


    jobject indexer = indexer_init();

    PG_TRY();
    {
        indexer_begin(indexer, indexPath);

        /* delete unnecessary data */

        Portal auditCursor = SPI_cursor_open_with_args(NULL, "select id from " AUDIT_TABLE " where not stored",
                0, NULL, NULL, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

        while(true)
        {
            SPI_cursor_fetch(auditCursor, true, SYNC_FETCH_SIZE);

            if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
                elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

            if(SPI_processed == 0)
                break;

            for(size_t i = 0; i < SPI_processed; i++)
            {
                HeapTuple tuple = SPI_tuptable->vals[i];

                Datum id = SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                indexer_delete(indexer, DatumGetInt32(id));
            }
        }

        SPI_cursor_close(auditCursor);


        /* convert new data */

        Portal compoundCursor = SPI_cursor_open_with_args(NULL, "select cmp.id, cmp.molfile from " COMPOUNDS_TABLE " cmp, "
                AUDIT_TABLE " aud where cmp.id = aud.id and aud.stored",
                0, NULL, NULL, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

        uint64 count = 0;

        while(true)
        {
            SPI_cursor_fetch(compoundCursor, true, SYNC_FETCH_SIZE);

            if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
                elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

            if(SPI_processed == 0)
                break;

            uint64 processed = SPI_processed;
            SPITupleTable *tuptable = SPI_tuptable;

            for(int i = 0; i < processed; i++)
            {
                HeapTuple tuple = tuptable->vals[i];


                Datum id = SPI_getbinval(tuple, tuptable->tupdesc, 1, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);


                Datum molfile = SPI_getbinval(tuple, tuptable->tupdesc, 2, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                VarChar *message = indexer_add(indexer, DatumGetInt32(id), DatumGetVarCharP(molfile));

                if(message != NULL)
                {
                    char *notice = text_to_cstring(message);
                    elog(NOTICE, "%i: %s", DatumGetInt32(id), notice);
                    pfree(notice);

                    if(SPI_execute_with_args("insert into " MOLECULE_ERRORS_TABLE " (compound, message) values ($1,$2)", 2,
                            (Oid[]) { INT4OID, TEXTOID },
                            (Datum[]) { id, PointerGetDatum(message) }, NULL, false, 0) != SPI_OK_INSERT)
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

        if(unlikely(SPI_exec("delete from " AUDIT_TABLE, 0) != SPI_OK_DELETE))
            elog(ERROR, "%s: SPI_exec() failed", __func__);


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


PG_FUNCTION_INFO_V1(sachem_cleanup);
Datum sachem_cleanup(PG_FUNCTION_ARGS)
{
    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    int indexNumber = -1;

    if(unlikely(SPI_exec("select id from " INDEX_TABLE, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(SPI_processed != 0)
    {
        if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
            elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

        bool isNullFlag;
        Datum number = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE) || isNullFlag)
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        indexNumber = DatumGetInt32(number);
    }

    SPI_finish();


    char *luceneIndexName = get_index_name(INDEX_PREFIX_NAME, indexNumber);


    int dirfd = -1;
    DIR *dp = NULL;

    if((dirfd = open(get_file_path(""), O_DIRECTORY)) == -1 && errno != ENOENT)
        elog(ERROR, "%s: open() failed", __func__);

    if(dirfd == -1)
        PG_RETURN_VOID();


    PG_TRY();
    {
        if((dp = fdopendir(dirfd)) == NULL)
            elog(ERROR, "%s: fdopendir() failed", __func__);


        struct dirent *ep;

        while((ep = readdir(dp)))
        {
            if(!strncmp(ep->d_name, INDEX_PREFIX_NAME, sizeof(INDEX_PREFIX_NAME) - 1))
            {
                if(!strcmp(ep->d_name, luceneIndexName))
                    continue;

                elog(NOTICE, "delete lucene index '%s'", ep->d_name);

                char *luceneIndexPath = get_file_path(ep->d_name);
                delete_directory(luceneIndexPath);
            }
            else if(strcmp(ep->d_name, ".") && strcmp(ep->d_name, ".."))
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
