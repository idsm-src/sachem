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
static jmethodID searchMethod;
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

        setFolderMethod = (*env)->GetMethodID(env, luceneClass, "setFolder", "(Ljava/lang/String;I)V");
        java_check_exception(__func__);

        beginMethod = (*env)->GetMethodID(env, luceneClass, "begin", "()V");
        java_check_exception(__func__);

        addMethod = (*env)->GetMethodID(env, luceneClass, "add", "(I[I)V");
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

        searchMethod = (*env)->GetMethodID(env, luceneClass, "search", "([I)[J");
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

    lucene->bitsetArray = NULL;
    lucene->bitsetWords = NULL;
}


void lucene_set_folder(Lucene *lucene, const char *path, int32_t maxId)
{
    jstring folder = NULL;

    PG_TRY();
    {
        folder = (*env)->NewStringUTF(env, path);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, lucene->instance, setFolderMethod, folder, maxId);
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


void lucene_add(Lucene *lucene, int32_t id, IntegerFingerprint fp)
{
    jintArray fpArray = NULL;

    PG_TRY();
    {
        fpArray = (jintArray) (*env)->NewIntArray(env, fp.size);
        java_check_exception(__func__);

        (*env)->SetIntArrayRegion(env, fpArray, 0, fp.size, (jint*) fp.data);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, lucene->instance, addMethod, (jint) id, fpArray);
        java_check_exception(__func__);

        JavaDeleteRef(fpArray);
    }
    PG_CATCH();
    {
        JavaDeleteRef(fpArray);

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


LuceneResultSet lucene_search(Lucene *lucene, IntegerFingerprint fp)
{
    jintArray fpArray = NULL;

    PG_TRY();
    {
        fpArray = (jintArray) (*env)->NewIntArray(env, fp.size);
        java_check_exception(__func__);

        (*env)->SetIntArrayRegion(env, fpArray, 0, fp.size, (jint*) fp.data);
        java_check_exception(__func__);

        lucene->bitsetArray = (jlongArray) (*env)->CallObjectMethod(env, lucene->instance, searchMethod, fpArray);
        java_check_exception(__func__);

        jsize size = (*env)->GetArrayLength(env, lucene->bitsetArray);
        lucene->bitsetWords = (*env)->GetLongArrayElements(env, lucene->bitsetArray, NULL);
        java_check_exception(__func__);

        bitset_init(&lucene->hits, (uint64_t *) lucene->bitsetWords, size);

        JavaDeleteRef(fpArray);
    }
    PG_CATCH();
    {
        JavaDeleteRef(fpArray);
        JavaDeleteLongArray(lucene->bitsetArray, lucene->bitsetWords, JNI_ABORT);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return (LuceneResultSet) { .possition = 0 };;
}


size_t lucene_get(Lucene *lucene, LuceneResultSet *resultSet, int32_t *buffer, size_t size)
{
    int ret = 0;

    while(size > 0)
    {
        int id = bitset_next_set_bit(&lucene->hits, resultSet->possition);
        resultSet->possition = id == -1 ? -1 : id + 1;

        if(id == -1)
        {
            JavaDeleteLongArray(lucene->bitsetArray, lucene->bitsetWords, JNI_ABORT);
            break;
        }

        *(buffer++) = id;

        ret++;
        size--;
    }

    return ret;
}


void lucene_fail(Lucene *lucene, LuceneResultSet *resultSet)
{
    JavaDeleteLongArray(lucene->bitsetArray, lucene->bitsetWords, JNI_ABORT);
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
