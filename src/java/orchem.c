#include <postgres.h>
#include <catalog/pg_type.h>
#include <stdbool.h>
#include "bitset.h"
#include "orchem.h"
#include "sachem.h"


static bool initialised = false;
static jclass byteArrayClass = NULL;

static jclass orchemSubstructureSearchClass = NULL;
static jclass orchemSubstructureQueryDataClass = NULL;
static jfieldID countsField = NULL;
static jfieldID fpField = NULL;
static jfieldID moleculeField = NULL;
static jfieldID restHField = NULL;
static jmethodID orchemSubstructureQueryDataMethod = NULL;

static jclass orchemSimilaritySearchClass = NULL;
static jmethodID orchemSimilarityQueryDataMethod = NULL;

static jclass orchemLoaderClass = NULL;
static jclass orchemLoaderDataClass = NULL;
static jfieldID orchemLoaderExceptionField = NULL;
static jfieldID orchemLoaderCountsField = NULL;
static jfieldID orchemLoaderFpField = NULL;
static jfieldID orchemLoaderMoleculeField = NULL;
static jmethodID orchemLoaderDataMethod = NULL;


void java_orchem_init(void)
{
    if(likely(initialised))
        return;

    java_init();


    byteArrayClass = (*env)->FindClass(env, "[B");
    java_check_exception(__func__);


    orchemSubstructureSearchClass = (*env)->FindClass(env, "cz/iocb/sachem/search/OrchemSubstructureSearch");
    java_check_exception(__func__);

    orchemSubstructureQueryDataClass = (*env)->FindClass(env, "cz/iocb/sachem/search/OrchemSubstructureSearch$OrchemQueryData");
    java_check_exception(__func__);

    countsField = (*env)->GetFieldID(env, orchemSubstructureQueryDataClass, "counts", "[S");
    java_check_exception(__func__);

    fpField = (*env)->GetFieldID(env, orchemSubstructureQueryDataClass, "fp", "[S");
    java_check_exception(__func__);

    moleculeField = (*env)->GetFieldID(env, orchemSubstructureQueryDataClass, "molecule", "[B");
    java_check_exception(__func__);

    restHField = (*env)->GetFieldID(env, orchemSubstructureQueryDataClass, "restH", "[Z");
    java_check_exception(__func__);

    orchemSubstructureQueryDataMethod = (*env)->GetStaticMethodID(env, orchemSubstructureSearchClass, "getQueryData", "([BIZZ)[Lcz/iocb/sachem/search/OrchemSubstructureSearch$OrchemQueryData;");
    java_check_exception(__func__);


    orchemSimilaritySearchClass = (*env)->FindClass(env, "cz/iocb/sachem/search/OrchemSimilaritySearch");
    java_check_exception(__func__);

    orchemSimilarityQueryDataMethod = (*env)->GetStaticMethodID(env, orchemSimilaritySearchClass, "getQueryData", "([BI)[J");
    java_check_exception(__func__);


    orchemLoaderClass = (*env)->FindClass(env, "cz/iocb/sachem/search/OrchemLoader");
    java_check_exception(__func__);

    orchemLoaderDataClass = (*env)->FindClass(env, "cz/iocb/sachem/search/OrchemLoader$OrchemData");
    java_check_exception(__func__);

    orchemLoaderExceptionField = (*env)->GetFieldID(env, orchemLoaderDataClass, "exception", "Ljava/lang/String;");
    java_check_exception(__func__);

    orchemLoaderCountsField = (*env)->GetFieldID(env, orchemLoaderDataClass, "counts", "[S");
    java_check_exception(__func__);

    orchemLoaderFpField = (*env)->GetFieldID(env, orchemLoaderDataClass, "fp", "[J");
    java_check_exception(__func__);

    orchemLoaderMoleculeField = (*env)->GetFieldID(env, orchemLoaderDataClass, "molecule", "[B");
    java_check_exception(__func__);

    orchemLoaderDataMethod = (*env)->GetStaticMethodID(env, orchemLoaderClass, "getIndexData", "([[B)[Lcz/iocb/sachem/search/OrchemLoader$OrchemData;");
    java_check_exception(__func__);


    initialised = true;
}


int java_orchem_parse_substructure_query(OrchemSubstructureQueryData **data, char* query, size_t queryLength, int32_t type, bool implicitHydrogens, bool tautomers)
{
    jbyteArray queryArg = NULL;
    jobjectArray result = NULL;
    jobject element = NULL;
    jshortArray countsArray = NULL;
    jshortArray fpArray = NULL;
    jbyteArray moleculeArray = NULL;
    jbooleanArray  restHArray = NULL;
    jshort *counts = NULL;
    jshort *fp = NULL;
    jbyte *molecule = NULL;
    jboolean *restH = NULL;
    jsize length = -1;


    PG_TRY();
    {
        queryArg = (*env)->NewByteArray(env, queryLength);
        java_check_exception(__func__);

        (*env)->SetByteArrayRegion(env, queryArg, 0, queryLength, (jbyte *) query);


        result = (jobjectArray) (*env)->CallStaticObjectMethod(env, orchemSubstructureSearchClass, orchemSubstructureQueryDataMethod, queryArg, (jint) type, (jboolean) implicitHydrogens, (jboolean) tautomers);
        java_check_exception(__func__);

        JavaDeleteRef(queryArg);


        length = (*env)->GetArrayLength(env, result);
        OrchemSubstructureQueryData *results = (OrchemSubstructureQueryData *) palloc(length * sizeof(OrchemSubstructureQueryData));

        for(int i = 0; i < length; i++)
        {
            element = (*env)->GetObjectArrayElement(env, result, i);

            countsArray = (jshortArray) (*env)->GetObjectField(env, element, countsField);
            fpArray = (jshortArray) (*env)->GetObjectField(env, element, fpField);
            moleculeArray = (jbyteArray)  (*env)->GetObjectField(env, element, moleculeField);
            restHArray = (jbooleanArray)   (*env)->GetObjectField(env, element, restHField);

            jsize countsSize = (*env)->GetArrayLength(env, countsArray);
            jsize fpSize = (*env)->GetArrayLength(env, fpArray);
            jsize moleculeSize = (*env)->GetArrayLength(env, moleculeArray);
            jsize restHSize = restHArray ? (*env)->GetArrayLength(env, restHArray) : -1;

            counts = (*env)->GetShortArrayElements(env, countsArray, NULL);
            java_check_exception(__func__);

            fp = (*env)->GetShortArrayElements(env, fpArray, NULL);
            java_check_exception(__func__);

            molecule = (*env)->GetByteArrayElements(env, moleculeArray, NULL);
            java_check_exception(__func__);

            restH = restHArray ? (*env)->GetBooleanArrayElements (env, restHArray, NULL) : 0;
            java_check_exception(__func__);


            results[i].counts = (jshort *) palloc(countsSize * sizeof(jshort));
            memcpy(results[i].counts, counts, countsSize * sizeof(jshort));

            results[i].fp = (jshort *) palloc(fpSize * sizeof(jshort));
            memcpy(results[i].fp, fp, fpSize * sizeof(jshort));

            results[i].molecule = (uint8_t *) palloc(moleculeSize);
            memcpy(results[i].molecule, molecule, moleculeSize);


            if(restHArray)
            {
                results[i].restH = (bool *) palloc(restHSize * sizeof(bool));
                memcpy(results[i].restH, restH, restHSize * sizeof(bool));
            }
            else
            {
                results[i].restH = NULL;
            }

            results[i].fpLength = fpSize;

            JavaDeleteShortArray(countsArray, counts, JNI_ABORT);
            JavaDeleteShortArray(fpArray, fp, JNI_ABORT);
            JavaDeleteByteArray(moleculeArray, molecule, JNI_ABORT);
            JavaDeleteBooleanArray(restHArray, restH, JNI_ABORT);

            JavaDeleteRef(element);
        }

        JavaDeleteRef(result);

        *data = results;
    }
    PG_CATCH();
    {
        JavaDeleteRef(queryArg);
        JavaDeleteRef(result);
        JavaDeleteRef(element);
        JavaDeleteShortArray(countsArray, counts, JNI_ABORT);
        JavaDeleteShortArray(fpArray, fp, JNI_ABORT);
        JavaDeleteByteArray(moleculeArray, molecule, JNI_ABORT);
        JavaDeleteBooleanArray(restHArray, restH, JNI_ABORT);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return length;
}


int java_orchem_parse_similarity_query(uint64_t **data, char* query, size_t queryLength, int32_t type)
{
    jbyteArray queryArg = NULL;
    jlongArray result = NULL;
    jlong *words = NULL;
    jsize length = -1;


    PG_TRY();
    {
        queryArg = (*env)->NewByteArray(env, queryLength);
        java_check_exception(__func__);

        (*env)->SetByteArrayRegion(env, queryArg, 0, queryLength, (jbyte *) query);

        result = (jlongArray) (*env)->CallStaticObjectMethod(env, orchemSimilaritySearchClass, orchemSimilarityQueryDataMethod, queryArg, (jint) type);
        java_check_exception(__func__);

        JavaDeleteRef(queryArg);

        if(result != NULL)
        {
            length = (*env)->GetArrayLength(env, result);
            words = (*env)->GetLongArrayElements(env, result, NULL);
            java_check_exception(__func__);

            *data = (uint64_t *) palloc(length * sizeof(jlong));
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
        JavaDeleteLongArray(result, words, JNI_ABORT);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return length;
}


void java_orchem_parse_data(size_t count, VarChar **molfiles, OrchemLoaderData *data)
{
    jbyteArray molfileArrayArg = NULL;
    jobjectArray resultArray = NULL;
    jbyteArray molfileArg = NULL;
    jobject resultElement = NULL;
    jstring exception = NULL;
    jshortArray countsArray = NULL;
    jlongArray fpArray = NULL;
    jbyteArray moleculeArray = NULL;
    jshort *counts = NULL;
    jlong *fp = NULL;
    jbyte *molecule = NULL;


    PG_TRY();
    {
        molfileArrayArg = (*env)->NewObjectArray(env, count, byteArrayClass, NULL);
        java_check_exception(__func__);

        for(size_t i = 0; i < count; i++)
        {
            int length = VARSIZE(molfiles[i]) - VARHDRSZ;
            molfileArg = (*env)->NewByteArray(env, length);
            java_check_exception(__func__);

            (*env)->SetByteArrayRegion(env, molfileArg, 0, length, (jbyte *) VARDATA(molfiles[i]));
            (*env)->SetObjectArrayElement(env, molfileArrayArg, i, molfileArg);

            JavaDeleteRef(molfileArg);
        }


        resultArray = (*env)->CallStaticObjectMethod(env, orchemLoaderClass, orchemLoaderDataMethod, molfileArrayArg);
        java_check_exception(__func__);

        JavaDeleteRef(molfileArrayArg);


        for(size_t i = 0; i < count; i++)
        {
            resultElement = (*env)->GetObjectArrayElement(env, resultArray, i);
            exception = (*env)->GetObjectField(env, resultElement, orchemLoaderExceptionField);

            if(exception)
            {
                jboolean isCopy;
                const char *message = (*env)->GetStringUTFChars(env, exception, &isCopy);
                data[i].error = cstring_to_text(message);
                (*env)->ReleaseStringUTFChars(env, exception, message);
                JavaDeleteRef(exception);

                data[i].bitCount = -1;
                data[i].fp = NULL;
                data[i].counts = NULL;
                data[i].molecule = NULL;
            }
            else
            {
                countsArray = (jshortArray) (*env)->GetObjectField(env, resultElement, orchemLoaderCountsField);
                fpArray = (jlongArray) (*env)->GetObjectField(env, resultElement, orchemLoaderFpField);
                moleculeArray = (jbyteArray)  (*env)->GetObjectField(env, resultElement, orchemLoaderMoleculeField);

                jsize countsSize = (*env)->GetArrayLength(env, countsArray);
                jsize fpSize = (*env)->GetArrayLength(env, fpArray);
                jsize moleculeSize = (*env)->GetArrayLength(env, moleculeArray);

                counts = (*env)->GetShortArrayElements(env, countsArray, NULL);
                java_check_exception(__func__);

                fp = (*env)->GetLongArrayElements(env, fpArray, NULL);
                java_check_exception(__func__);

                molecule = (*env)->GetByteArrayElements(env, moleculeArray, NULL);
                java_check_exception(__func__);


                data[i].fp = (ArrayType *) palloc(fpSize * sizeof(uint64_t) + ARR_OVERHEAD_NONULLS(1));
                data[i].fp->ndim = 1;
                data[i].fp->dataoffset = 0;
                data[i].fp->elemtype = INT8OID;
                memcpy(ARR_DATA_PTR(data[i].fp), fp, fpSize * sizeof(uint64_t));
                *(ARR_DIMS(data[i].fp)) = fpSize;
                *(ARR_LBOUND(data[i].fp)) = 1;
                SET_VARSIZE(data[i].fp, fpSize * sizeof(uint64_t) + ARR_OVERHEAD_NONULLS(1));

                BitSet bitset;
                bitset_init(&bitset, (uint64_t *) fp, fpSize);
                data[i].bitCount = bitset_cardinality(&bitset);

                data[i].counts = (ArrayType *) palloc(countsSize * sizeof(int16) + ARR_OVERHEAD_NONULLS(1));
                data[i].counts->ndim = 1;
                data[i].counts->dataoffset = 0;
                data[i].counts->elemtype = INT2OID;
                memcpy(ARR_DATA_PTR(data[i].counts), counts, countsSize * sizeof(int16));
                *(ARR_DIMS(data[i].counts)) = countsSize;
                *(ARR_LBOUND(data[i].counts)) = 1;
                SET_VARSIZE(data[i].counts, countsSize * sizeof(int16) + ARR_OVERHEAD_NONULLS(1));

                data[i].molecule = (bytea *) palloc(VARHDRSZ + moleculeSize);
                SET_VARSIZE(data[i].molecule, VARHDRSZ + moleculeSize);
                memcpy(VARDATA(data[i].molecule), molecule, moleculeSize);

                data[i].error = NULL;


                JavaDeleteShortArray(countsArray, counts, JNI_ABORT);
                JavaDeleteLongArray(fpArray, fp, JNI_ABORT);
                JavaDeleteByteArray(moleculeArray, molecule, JNI_ABORT);
            }

            JavaDeleteRef(resultElement);
        }

        JavaDeleteRef(resultArray);
    }
    PG_CATCH();
    {
        JavaDeleteRef(molfileArrayArg);
        JavaDeleteRef(resultArray);
        JavaDeleteRef(molfileArg);
        JavaDeleteRef(resultElement);
        JavaDeleteRef(exception);
        JavaDeleteShortArray(countsArray, counts, JNI_ABORT);
        JavaDeleteLongArray(fpArray, fp, JNI_ABORT);
        JavaDeleteByteArray(moleculeArray, molecule, JNI_ABORT);

        PG_RE_THROW();
    }
    PG_END_TRY();
}
