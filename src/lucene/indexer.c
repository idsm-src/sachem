#include <postgres.h>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "bitset.h"
#include "sachem.h"
#include "indexer.h"
#include "java/java.h"


static bool indexerInitialized = false;
static jclass indexerClass;
static jmethodID constructor;
static jmethodID beginMethod;
static jmethodID addMethod;
static jmethodID addIndexMethod;
static jmethodID deleteMethod;
static jmethodID optimizeMethod;
static jmethodID commitMethod;
static jmethodID rollbackMethod;


static void lucene_indexer_java_init()
{
    if(likely(indexerInitialized))
        return;

    java_init();

    indexerClass = (*env)->FindClass(env, "cz/iocb/sachem/lucene/Indexer");
    java_check_exception(__func__);

    constructor = (*env)->GetMethodID(env, indexerClass, "<init>", "()V");
    java_check_exception(__func__);

    beginMethod = (*env)->GetMethodID(env, indexerClass, "begin", "(Ljava/lang/String;)V");
    java_check_exception(__func__);

    addMethod = (*env)->GetMethodID(env, indexerClass, "add", "(I[I[I)V");
    java_check_exception(__func__);

    addIndexMethod = (*env)->GetMethodID(env, indexerClass, "addIndex", "(Ljava/lang/String;)V");
    java_check_exception(__func__);

    deleteMethod = (*env)->GetMethodID(env, indexerClass, "delete", "(I)V");
    java_check_exception(__func__);

    optimizeMethod = (*env)->GetMethodID(env, indexerClass, "optimize", "()V");
    java_check_exception(__func__);

    commitMethod = (*env)->GetMethodID(env, indexerClass, "commit", "()V");
    java_check_exception(__func__);

    rollbackMethod = (*env)->GetMethodID(env, indexerClass, "rollback", "()V");
    java_check_exception(__func__);

    indexerInitialized = true;
}


void lucene_indexer_init(LuceneIndexer *lucene)
{
    lucene_indexer_java_init();

    lucene->instance = (*env)->NewObject(env, indexerClass, constructor);
    java_check_exception(__func__);
}


void lucene_indexer_terminate(LuceneIndexer *lucene)
{
    JavaDeleteRef(lucene->instance);
}


void lucene_indexer_begin(LuceneIndexer *lucene, const char *path)
{
    jstring folder = NULL;

    PG_TRY();
    {
        folder = (*env)->NewStringUTF(env, path);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, lucene->instance, beginMethod, folder);
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


void lucene_indexer_add(LuceneIndexer *lucene, int32_t id, IntegerFingerprint subfp, IntegerFingerprint simfp)
{
    jintArray subfpArray = NULL;
    jintArray simfpArray = NULL;

    PG_TRY();
    {
        subfpArray = (jintArray) (*env)->NewIntArray(env, subfp.size);
        java_check_exception(__func__);

        (*env)->SetIntArrayRegion(env, subfpArray, 0, subfp.size, (jint*) subfp.data);
        java_check_exception(__func__);


        simfpArray = (jintArray) (*env)->NewIntArray(env, simfp.size);
        java_check_exception(__func__);

        (*env)->SetIntArrayRegion(env, simfpArray, 0, simfp.size, (jint*) simfp.data);
        java_check_exception(__func__);


        (*env)->CallVoidMethod(env, lucene->instance, addMethod, (jint) id, subfpArray, simfpArray);
        java_check_exception(__func__);

        JavaDeleteRef(subfpArray);
        JavaDeleteRef(simfpArray);
    }
    PG_CATCH();
    {
        JavaDeleteRef(subfpArray);
        JavaDeleteRef(simfpArray);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


void lucene_indexer_add_index(LuceneIndexer *lucene, const char *path)
{
    jstring folder = NULL;

    PG_TRY();
    {
        folder = (*env)->NewStringUTF(env, path);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, lucene->instance, addIndexMethod, folder);
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


void lucene_indexer_delete(LuceneIndexer *lucene, int32_t id)
{
    (*env)->CallVoidMethod(env, lucene->instance, deleteMethod, (jint) id);
    java_check_exception(__func__);
}


void lucene_indexer_optimize(LuceneIndexer *lucene)
{
    (*env)->CallVoidMethod(env, lucene->instance, optimizeMethod);
    java_check_exception(__func__);
}


void lucene_indexer_commit(LuceneIndexer *lucene)
{
    (*env)->CallVoidMethod(env, lucene->instance, commitMethod);
    java_check_exception(__func__);
}


void lucene_indexer_rollback(LuceneIndexer *lucene)
{
    (*env)->CallVoidMethod(env, lucene->instance, rollbackMethod);
    java_check_exception(__func__);
}


void lucene_indexer_link_directory(const char *oldPath, const char *newPath)
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


void lucene_indexer_delete_directory(const char *path)
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
