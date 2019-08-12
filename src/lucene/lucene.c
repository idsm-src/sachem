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
static jmethodID subsearchMethod;
static jmethodID simsearchMethod;
static jclass scoreHitClass;
static jfieldID idField;
static jfieldID scoreField;


void lucene_java_init()
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
}


void lucene_init(Lucene *lucene)
{
    lucene_java_init();

    lucene->instance = (*env)->NewObject(env, luceneClass, constructor);
    java_check_exception(__func__);
}


void lucene_terminate(Lucene *lucene)
{
    JavaDeleteRef(lucene->instance);
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
