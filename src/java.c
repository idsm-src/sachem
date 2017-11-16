#include <jni.h>
#include <postgres.h>
#include <libpq/pqsignal.h>
#include <storage/ipc.h>
#include <tcop/tcopprot.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include "java.h"
#include "sachem.h"


#define JavaDeleteRef(ref) \
    do \
    { \
        if(likely(ref != NULL)) \
        (*env)->DeleteLocalRef(env, ref); \
        ref = NULL; \
    } while(0)


#define JavaDeleteBooleanArray(array, buffer, mode) \
    do { \
        if(likely(array != NULL && buffer != NULL)) (*env)->ReleaseBooleanArrayElements(env, array, buffer, mode); \
        if(likely(array != NULL)) (*env)->DeleteLocalRef(env, array); \
        array = NULL; \
        buffer = NULL; \
    } while(0)


#define JavaDeleteByteArray(array, buffer, mode) \
    do { \
        if(likely(array != NULL && buffer != NULL)) (*env)->ReleaseByteArrayElements(env, array, buffer, mode); \
        if(likely(array != NULL)) (*env)->DeleteLocalRef(env, array); \
        array = NULL; \
        buffer = NULL; \
    } while(0)


#define JavaDeleteShortArray(array, buffer, mode) \
    do { \
        if(likely(array != NULL && buffer != NULL)) (*env)->ReleaseShortArrayElements(env, array, buffer, mode); \
        if(likely(array != NULL)) (*env)->DeleteLocalRef(env, array); \
        array = NULL; \
        buffer = NULL; \
    } while(0)


#define JavaDeleteLongArray(array, buffer, mode) \
    do { \
        if(likely(array != NULL && buffer != NULL)) (*env)->ReleaseLongArrayElements(env, array, buffer, mode); \
        if(likely(array != NULL)) (*env)->DeleteLocalRef(env, array); \
        array = NULL; \
        buffer = NULL; \
    } while(0)


static bool initialised = false;
static JavaVM* jvm = NULL;
static JNIEnv* env = NULL;
static jclass substructureSearchClass = NULL;
static jclass substructureQueryDataClass = NULL;
static jmethodID queryDataMethod = NULL;
static jfieldID countsField = NULL;
static jfieldID fpField = NULL;
static jfieldID atomsField = NULL;
static jfieldID bondsField = NULL;
static jfieldID restHField = NULL;
static jclass exceptionClass = NULL;
static jmethodID toStringMethod = NULL;
static jclass similaritySearchClass = NULL;
static jmethodID similarityQueryDataMethod = NULL;
static pthread_t mainThread;


static inline void java_check_exception(const char *str)
{
    jthrowable exception = (*env)->ExceptionOccurred(env);

    if(exception == NULL)
        return;

    (*env)->ExceptionDescribe(env);
    (*env)->ExceptionClear(env);

    jstring message = (jstring)(*env)->CallObjectMethod(env, exception, toStringMethod);

    const char *mstr = (*env)->GetStringUTFChars(env, message, NULL);
    char *error = palloc(strlen(mstr) + 1);
    strcpy(error, mstr);

    (*env)->ReleaseStringUTFChars(env, message, mstr);
    JavaDeleteRef(message);
    JavaDeleteRef(exception);

    elog(ERROR, "%s: java error: %s", str, error);
}


static void java_statement_cancel_handler(int signum)
{
    if(pthread_self() == mainThread)
        StatementCancelHandler(signum);
    else
        pthread_kill(mainThread, signum);
}


static void java_die_handler(int signum)
{
    if(pthread_self() == mainThread)
        die(signum);
    else
        pthread_kill(mainThread, signum);
}


static void java_quick_die_handler(int signum)
{
    if(pthread_self() == mainThread)
        quickdie(signum);
    else
        pthread_kill(mainThread, signum);
}


void java_module_init(void)
{
    JavaVMInitArgs args = (JavaVMInitArgs) {
        .version = JNI_VERSION_1_8,
        .nOptions = 2,
        .options = (JavaVMOption []) {
            (JavaVMOption) { .optionString = "-Djava.class.path=" JARDIR "/orchem.jar:" JARDIR "/cdk-2.0.jar"},
            (JavaVMOption) { .optionString = "-Dinchi.path=" BINDIR "/inchi-1"}
        },
        .ignoreUnrecognized = JNI_FALSE
    };

    pqsignal(SIGINT, SIG_IGN);
    pqsignal(SIGTERM, SIG_IGN);
    pqsignal(SIGQUIT, SIG_IGN);

    JNI_CreateJavaVM(&(jvm), (void **) &(env), &args);

    mainThread = pthread_self();
    pqsignal(SIGINT,  java_statement_cancel_handler);
    pqsignal(SIGTERM, java_die_handler);
    pqsignal(SIGQUIT, java_quick_die_handler);

    if(jvm == NULL || env == NULL)
        elog(ERROR, "cannot initialize JVM");


    substructureSearchClass = (*env)->FindClass(env, "cz/iocb/orchem/search/SubstructureSearch");
    java_check_exception("java_module_init()");

    substructureQueryDataClass = (*env)->FindClass(env, "cz/iocb/orchem/search/SubstructureSearch$SubstructureQueryData");
    java_check_exception("java_module_init()");

    queryDataMethod = (*env)->GetStaticMethodID(env, substructureSearchClass, "getQueryData", "([BLjava/lang/String;Z)[Lcz/iocb/orchem/search/SubstructureSearch$SubstructureQueryData;");
    java_check_exception("java_module_init()");

    countsField = (*env)->GetFieldID(env, substructureQueryDataClass, "counts", "[S");
    java_check_exception("java_module_init()");

    fpField = (*env)->GetFieldID(env, substructureQueryDataClass, "fp", "[S");
    java_check_exception("java_module_init()");

    atomsField = (*env)->GetFieldID(env, substructureQueryDataClass, "atoms", "[B");
    java_check_exception("java_module_init()");

    bondsField = (*env)->GetFieldID(env, substructureQueryDataClass, "bonds", "[B");
    java_check_exception("java_module_init()");

    restHField = (*env)->GetFieldID(env, substructureQueryDataClass, "restH", "[Z");
    java_check_exception("java_module_init()");

    exceptionClass = (*env)->FindClass(env, "java/lang/Throwable");
    java_check_exception("java_module_init()");

    toStringMethod = (*env)->GetMethodID(env, exceptionClass, "toString", "()Ljava/lang/String;");
    java_check_exception("java_module_init()");


    similaritySearchClass = (*env)->FindClass(env, "cz/iocb/orchem/search/SimilaritySearch");
    java_check_exception("java_module_init()");

    similarityQueryDataMethod = (*env)->GetStaticMethodID(env, similaritySearchClass, "getQueryData", "([BLjava/lang/String;)[J");
    java_check_exception("java_module_init()");


    initialised = true;
}


void java_module_finish(void)
{
    if(jvm != NULL)
        (*jvm)->DestroyJavaVM(jvm);
}


int java_parse_query(QueryData **data, char* query, size_t queryLength, char *type, bool tautomers)
{
    if(initialised == false)
        elog(ERROR, "java module is not properly initialized");


    jbyteArray queryArg = NULL;
    jstring typeArg = NULL;
    jobjectArray result = NULL;
    jobject element = NULL;
    jshortArray countsArray = NULL;
    jshortArray fpArray = NULL;
    jbyteArray atomsArray = NULL;
    jbyteArray bondsArray = NULL;
    jbooleanArray  restHArray = NULL;
    jshort *counts = NULL;
    jshort *fp = NULL;
    jbyte *atoms = NULL;
    jbyte *bonds = NULL;
    jboolean *restH = NULL;
    jsize length = -1;


    PG_TRY();
    {
        queryArg = (*env)->NewByteArray(env, queryLength);
        java_check_exception("java_parse_query()");

        (*env)->SetByteArrayRegion(env, queryArg, 0, queryLength, (jbyte*) query);


        typeArg = (*env)->NewStringUTF(env, type);
        java_check_exception("java_parse_query()");


        result = (jobjectArray) (*env)->CallStaticObjectMethod(env, substructureSearchClass, queryDataMethod, queryArg, typeArg, (jboolean) tautomers);
        java_check_exception("java_parse_query()");

        JavaDeleteRef(queryArg);
        JavaDeleteRef(typeArg);


        length = (*env)->GetArrayLength(env, result);
        QueryData *results = palloc(length * sizeof(QueryData));

        for(int i = 0; i < length; i++)
        {
            element = (*env)->GetObjectArrayElement(env, result, i);

            countsArray = (jshortArray) (*env)->GetObjectField(env, element, countsField);
            fpArray = (jshortArray) (*env)->GetObjectField(env, element, fpField);
            atomsArray = (jbyteArray)  (*env)->GetObjectField(env, element, atomsField);
            bondsArray = (jbyteArray)  (*env)->GetObjectField(env, element, bondsField);
            restHArray = (jbooleanArray)   (*env)->GetObjectField(env, element, restHField);

            jsize countsSize = (*env)->GetArrayLength(env, countsArray);
            jsize fpSize = (*env)->GetArrayLength(env, fpArray);
            jsize atomsSize = (*env)->GetArrayLength(env, atomsArray);
            jsize bondsSize = (*env)->GetArrayLength(env, bondsArray);
            jsize restHSize = restHArray ? (*env)->GetArrayLength(env, restHArray) : -1;

            counts = (*env)->GetShortArrayElements(env, countsArray, NULL);
            java_check_exception("java_parse_query()");

            fp = (*env)->GetShortArrayElements(env, fpArray, NULL);
            java_check_exception("java_parse_query()");

            atoms = (*env)->GetByteArrayElements(env, atomsArray, NULL);
            java_check_exception("java_parse_query()");

            bonds = (*env)->GetByteArrayElements(env, bondsArray, NULL);
            java_check_exception("java_parse_query()");

            restH = restHArray ? (*env)->GetBooleanArrayElements (env, restHArray, NULL) : 0;
            java_check_exception("java_parse_query()");


            results[i].counts = palloc(countsSize * sizeof(jshort));
            memcpy(results[i].counts, counts, countsSize * sizeof(jshort));

            results[i].fp = palloc(fpSize * sizeof(jshort));
            memcpy(results[i].fp, fp, fpSize * sizeof(jshort));

            results[i].atoms = palloc(atomsSize);
            memcpy(results[i].atoms, atoms, atomsSize);

            results[i].bonds = palloc(bondsSize);
            memcpy(results[i].bonds, bonds, bondsSize);


            if(restHArray)
            {
                results[i].restH = palloc(restHSize * sizeof(bool));
                memcpy(results[i].restH, restH, restHSize * sizeof(bool));
            }
            else
            {
                results[i].restH = NULL;
            }

            results[i].fpLength = fpSize;
            results[i].atomLength = atomsSize;
            results[i].bondLength = bondsSize;

            JavaDeleteShortArray(countsArray, counts, JNI_ABORT);
            JavaDeleteShortArray(fpArray, fp, JNI_ABORT);
            JavaDeleteByteArray(atomsArray, atoms, JNI_ABORT);
            JavaDeleteByteArray(bondsArray, bonds, JNI_ABORT);
            JavaDeleteBooleanArray(restHArray, restH, JNI_ABORT);

            JavaDeleteRef(element);
        }

        JavaDeleteRef(result);

        *data = results;
    }
    PG_CATCH();
    {
        JavaDeleteRef(queryArg);
        JavaDeleteRef(typeArg);
        JavaDeleteRef(result);
        JavaDeleteRef(element);
        JavaDeleteShortArray(countsArray, counts, JNI_ABORT);
        JavaDeleteShortArray(fpArray, fp, JNI_ABORT);
        JavaDeleteByteArray(atomsArray, atoms, JNI_ABORT);
        JavaDeleteByteArray(bondsArray, bonds, JNI_ABORT);
        JavaDeleteBooleanArray(restHArray, restH, JNI_ABORT);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return length;
}


int java_parse_similarity_query(uint64_t **data, char* query, size_t queryLength, char *type)
{
    if(initialised == false)
        elog(ERROR, "java module is not properly initialized");

    jbyteArray queryArg = NULL;
    jstring typeArg = NULL;
    jlongArray result = NULL;
    jlong *words = NULL;
    jsize length = -1;


    PG_TRY();
    {
        queryArg = (*env)->NewByteArray(env, queryLength);
        java_check_exception("java_parse_similarity_query()");

        (*env)->SetByteArrayRegion(env, queryArg, 0, queryLength, (jbyte*) query);

        typeArg = (*env)->NewStringUTF(env, type);
        java_check_exception("java_parse_similarity_query()");

        result = (jlongArray) (*env)->CallStaticObjectMethod(env, similaritySearchClass, similarityQueryDataMethod, queryArg, typeArg);
        java_check_exception("java_parse_similarity_query()");

        JavaDeleteRef(queryArg);
        JavaDeleteRef(typeArg);

        if(result != NULL)
        {
            length = (*env)->GetArrayLength(env, result);
            words = (*env)->GetLongArrayElements(env, result, NULL);
            java_check_exception("java_parse_similarity_query()");

            *data = palloc(length * sizeof(jlong));
            memcpy(*data, words, length * sizeof(jlong));

            JavaDeleteLongArray(result, words, JNI_ABORT);
        }
        else
        {
            *data = NULL;
        }
    }
    PG_CATCH();
    {
        JavaDeleteRef(queryArg);
        JavaDeleteRef(typeArg);
        JavaDeleteLongArray(result, words, JNI_ABORT);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return length;
}
