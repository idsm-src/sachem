#include <postgres.h>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "bitset.h"
#include "sachem.h"
#include "lucene.h"
#include "java/java.h"


static bool luceneInitialized = false;
static jclass luceneClass;
static jmethodID constructor;
static jmethodID setFolderMethod;
static jmethodID beginMethod;
static jmethodID addMethod;
static jmethodID addIndexMethod;
static jmethodID deleteMethod;
static jmethodID optimizeMethod;
static jmethodID commitMethod;
static jmethodID rollbackMethod;
static jmethodID subsearchMethod;
static jmethodID simsearchMethod;
static jclass scoreHitClass;
static jfieldID idField;
static jfieldID scoreField;


void lucene_init(Lucene *lucene)
{
    if(luceneInitialized == false)
    {
        java_init();

        luceneClass = (*env)->FindClass(env, "cz/iocb/sachem/lucene/Lucene");
        java_check_exception(__func__);

        constructor = (*env)->GetMethodID(env, luceneClass, "<init>", "()V");
        java_check_exception(__func__);

        setFolderMethod = (*env)->GetMethodID(env, luceneClass, "setFolder", "(Ljava/lang/String;)V");
        java_check_exception(__func__);

        beginMethod = (*env)->GetMethodID(env, luceneClass, "begin", "()V");
        java_check_exception(__func__);

        addMethod = (*env)->GetMethodID(env, luceneClass, "add", "(I[I[I)V");
        java_check_exception(__func__);

        addIndexMethod = (*env)->GetMethodID(env, luceneClass, "addIndex", "(Ljava/lang/String;)V");
        java_check_exception(__func__);

        deleteMethod = (*env)->GetMethodID(env, luceneClass, "delete", "(I)V");
        java_check_exception(__func__);

        optimizeMethod = (*env)->GetMethodID(env, luceneClass, "optimize", "()V");
        java_check_exception(__func__);

        commitMethod = (*env)->GetMethodID(env, luceneClass, "commit", "()V");
        java_check_exception(__func__);

        rollbackMethod = (*env)->GetMethodID(env, luceneClass, "rollback", "()V");
        java_check_exception(__func__);

        subsearchMethod = (*env)->GetMethodID(env, luceneClass, "subsearch", "([II)[J");
        java_check_exception(__func__);

        simsearchMethod = (*env)->GetMethodID(env, luceneClass, "simsearch", "([IIF)[Lcz/iocb/sachem/lucene/ScoreHit;");
        java_check_exception(__func__);

        scoreHitClass = (*env)->FindClass(env, "cz/iocb/sachem/lucene/ScoreHit");
        java_check_exception(__func__);

        idField = (*env)->GetFieldID(env, scoreHitClass, "id", "I");
        java_check_exception(__func__);

        scoreField = (*env)->GetFieldID(env, scoreHitClass, "score", "F");
        java_check_exception(__func__);

        luceneInitialized = true;
    }

    lucene->instance = (*env)->NewObject(env, luceneClass, constructor);
    java_check_exception(__func__);
}


void lucene_set_folder(Lucene *lucene, const char *path)
{
    jstring folder = NULL;

    PG_TRY();
    {
        folder = (*env)->NewStringUTF(env, path);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, lucene->instance, setFolderMethod, folder);
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


void lucene_begin(Lucene *lucene)
{
    (*env)->CallVoidMethod(env, lucene->instance, beginMethod);
    java_check_exception(__func__);
}


void lucene_add(Lucene *lucene, int32_t id, IntegerFingerprint subfp, IntegerFingerprint simfp)
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


void lucene_add_index(Lucene *lucene, const char *path)
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


void lucene_delete(Lucene *lucene, int32_t id)
{
    (*env)->CallVoidMethod(env, lucene->instance, deleteMethod, (jint) id);
    java_check_exception(__func__);
}


void lucene_optimize(Lucene *lucene)
{
    (*env)->CallVoidMethod(env, lucene->instance, optimizeMethod);
    java_check_exception(__func__);
}


void lucene_commit(Lucene *lucene)
{
    (*env)->CallVoidMethod(env, lucene->instance, commitMethod);
    java_check_exception(__func__);
}


void lucene_rollback(Lucene *lucene)
{
    (*env)->CallVoidMethod(env, lucene->instance, rollbackMethod);
    java_check_exception(__func__);
}


LuceneSubsearchResult lucene_subsearch_submit(Lucene *lucene, IntegerFingerprint fp, int32_t maxId)
{
    LuceneSubsearchResult result = { .bitsetArray = NULL, .bitsetWords = NULL };
    jintArray fpArray = NULL;

    PG_TRY();
    {
        fpArray = (jintArray) (*env)->NewIntArray(env, fp.size);
        java_check_exception(__func__);

        (*env)->SetIntArrayRegion(env, fpArray, 0, fp.size, (jint*) fp.data);
        java_check_exception(__func__);

        result.bitsetArray = (jlongArray) (*env)->CallObjectMethod(env, lucene->instance, subsearchMethod, fpArray, maxId);
        java_check_exception(__func__);

        jsize size = (*env)->GetArrayLength(env, result.bitsetArray);
        result.bitsetWords = (*env)->GetLongArrayElements(env, result.bitsetArray, NULL);
        java_check_exception(__func__);

        bitset_init(&result.hits, (uint64_t *) result.bitsetWords, size);

        JavaDeleteRef(fpArray);
    }
    PG_CATCH();
    {
        JavaDeleteRef(fpArray);
        JavaDeleteLongArray(result.bitsetArray, result.bitsetWords, JNI_ABORT);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return result;
}


size_t lucene_subsearch_get(Lucene *lucene, LuceneSubsearchResult *result, int32_t *buffer, size_t size)
{
    int ret = 0;

    while(size > 0)
    {
        int id = bitset_next_set_bit(&result->hits, result->possition);
        result->possition = id == -1 ? -1 : id + 1;

        if(id == -1)
        {
            JavaDeleteLongArray(result->bitsetArray, result->bitsetWords, JNI_ABORT);
            break;
        }

        *(buffer++) = id;

        ret++;
        size--;
    }

    return ret;
}


void lucene_subsearch_fail(Lucene *lucene, LuceneSubsearchResult *result)
{
    JavaDeleteLongArray(result->bitsetArray, result->bitsetWords, JNI_ABORT);
}


LuceneSimsearchResult lucene_simsearch_submit(Lucene *lucene, IntegerFingerprint fp, int32_t topN, float cutoff)
{
    jintArray fpArray = NULL;
    jobjectArray result = NULL;
    jsize count;

    PG_TRY();
    {
        fpArray = (jintArray) (*env)->NewIntArray(env, fp.size);
        java_check_exception(__func__);

        (*env)->SetIntArrayRegion(env, fpArray, 0, fp.size, (jint*) fp.data);
        java_check_exception(__func__);

        result = (jobjectArray) (*env)->CallObjectMethod(env, lucene->instance, simsearchMethod, fpArray, topN, cutoff);
        java_check_exception(__func__);

        count = (*env)->GetArrayLength(env, result);

        JavaDeleteRef(fpArray);
    }
    PG_CATCH();
    {
        JavaDeleteRef(fpArray);
        JavaDeleteRef(result);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return (LuceneSimsearchResult) { .count = count, .possition = 0, .result = result };
}


bool lucene_simsearch_get(Lucene *lucene, LuceneSimsearchResult *result, int32_t *id, float *score)
{
    if(result->possition == result->count)
    {
        JavaDeleteRef(result->result);
        return false;
    }


    jobject element = NULL;

    PG_TRY();
    {
        element = (*env)->GetObjectArrayElement(env, result->result, result->possition);
        *id = (*env)->GetIntField(env, element, idField);
        *score = (*env)->GetFloatField(env, element, scoreField);

        java_check_exception(__func__);

        JavaDeleteRef(element);
    }
    PG_CATCH();
    {
        JavaDeleteRef(element);

        PG_RE_THROW();
    }
    PG_END_TRY();

    result->possition++;
    return true;
}


void lucene_simsearch_fail(Lucene *lucene, LuceneSimsearchResult *result)
{
    JavaDeleteRef(result->result);
}


void lucene_link_directory(const char *oldPath, const char *newPath)
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


void lucene_delete_directory(const char *path)
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
