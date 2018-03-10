#include <postgres.h>
#include <catalog/pg_type.h>
#include <stdbool.h>
#include "sachem.h"
#include "parse.h"


static bool initialised = false;
static jclass byteArrayClass = NULL;

static jclass substructureSearchClass = NULL;
static jclass substructureQueryDataClass = NULL;
static jfieldID moleculeField = NULL;
static jfieldID restHField = NULL;
static jmethodID substructureQueryDataMethod = NULL;

static jclass lucyLoaderClass = NULL;
static jclass lucyLoaderDataClass = NULL;
static jfieldID lucyLoaderExceptionField = NULL;
static jfieldID lucyLoaderMoleculeField = NULL;
static jmethodID lucyLoaderDataMethod = NULL;


void java_parse_init(void)
{
    if(likely(initialised))
        return;

    java_init();


    byteArrayClass = (*env)->FindClass(env, "[B");
    java_check_exception(__func__);


    substructureSearchClass = (*env)->FindClass(env, "cz/iocb/sachem/search/SubstructureSearch");
    java_check_exception(__func__);

    substructureQueryDataClass = (*env)->FindClass(env, "cz/iocb/sachem/search/SubstructureSearch$QueryData");
    java_check_exception(__func__);

    moleculeField = (*env)->GetFieldID(env, substructureQueryDataClass, "molecule", "[B");
    java_check_exception(__func__);

    restHField = (*env)->GetFieldID(env, substructureQueryDataClass, "restH", "[Z");
    java_check_exception(__func__);

    substructureQueryDataMethod = (*env)->GetStaticMethodID(env, substructureSearchClass, "getQueryData", "([BIZZ)[Lcz/iocb/sachem/search/SubstructureSearch$QueryData;");
    java_check_exception(__func__);


    lucyLoaderClass = (*env)->FindClass(env, "cz/iocb/sachem/search/LucyLoader");
    java_check_exception(__func__);

    lucyLoaderDataClass = (*env)->FindClass(env, "cz/iocb/sachem/search/LucyLoader$LucyData");
    java_check_exception(__func__);

    lucyLoaderExceptionField = (*env)->GetFieldID(env, lucyLoaderDataClass, "exception", "Ljava/lang/String;");
    java_check_exception(__func__);

    lucyLoaderMoleculeField = (*env)->GetFieldID(env, lucyLoaderDataClass, "molecule", "[B");
    java_check_exception(__func__);

    lucyLoaderDataMethod = (*env)->GetStaticMethodID(env, lucyLoaderClass, "getIndexData", "([[B)[Lcz/iocb/sachem/search/LucyLoader$LucyData;");
    java_check_exception(__func__);


    initialised = true;
}


int java_parse_substructure_query(SubstructureQueryData **data, char* query, size_t queryLength, int32_t type, bool implicitHydrogens, bool tautomers)
{
    jbyteArray queryArg = NULL;
    jobjectArray result = NULL;
    jobject element = NULL;
    jbyteArray moleculeArray = NULL;
    jbooleanArray restHArray = NULL;
    jbyte *molecule = NULL;
    jboolean *restH = NULL;
    jsize length = -1;


    PG_TRY();
    {
        queryArg = (*env)->NewByteArray(env, queryLength);
        java_check_exception(__func__);

        (*env)->SetByteArrayRegion(env, queryArg, 0, queryLength, (jbyte *) query);


        result = (jobjectArray) (*env)->CallStaticObjectMethod(env, substructureSearchClass, substructureQueryDataMethod, queryArg, (jint) type, (jboolean) implicitHydrogens, (jboolean) tautomers);
        java_check_exception(__func__);

        JavaDeleteRef(queryArg);


        length = (*env)->GetArrayLength(env, result);
        SubstructureQueryData *results = (SubstructureQueryData *) palloc(length * sizeof(SubstructureQueryData));

        for(int i = 0; i < length; i++)
        {
            element = (*env)->GetObjectArrayElement(env, result, i);

            moleculeArray = (jbyteArray) (*env)->GetObjectField(env, element, moleculeField);
            restHArray = (jbooleanArray) (*env)->GetObjectField(env, element, restHField);

            jsize moleculeSize = (*env)->GetArrayLength(env, moleculeArray);
            jsize restHSize = restHArray ? (*env)->GetArrayLength(env, restHArray) : -1;

            molecule = (*env)->GetByteArrayElements(env, moleculeArray, NULL);
            java_check_exception(__func__);

            restH = restHArray ? (*env)->GetBooleanArrayElements (env, restHArray, NULL) : 0;
            java_check_exception(__func__);


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
        JavaDeleteByteArray(moleculeArray, molecule, JNI_ABORT);
        JavaDeleteBooleanArray(restHArray, restH, JNI_ABORT);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return length;
}


void java_parse_data(size_t count, VarChar **molfiles, LoaderData *data)
{
    jbyteArray molfileArrayArg = NULL;
    jobjectArray resultArray = NULL;
    jbyteArray molfileArg = NULL;
    jobject resultElement = NULL;
    jstring exception = NULL;
    jbyteArray moleculeArray = NULL;
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


        resultArray = (*env)->CallStaticObjectMethod(env, lucyLoaderClass, lucyLoaderDataMethod, molfileArrayArg);
        java_check_exception(__func__);

        JavaDeleteRef(molfileArrayArg);


        for(size_t i = 0; i < count; i++)
        {
            resultElement = (*env)->GetObjectArrayElement(env, resultArray, i);
            exception = (*env)->GetObjectField(env, resultElement, lucyLoaderExceptionField);

            if(exception)
            {
                jboolean isCopy;
                const char *message = (*env)->GetStringUTFChars(env, exception, &isCopy);
                data[i].error = cstring_to_text(message);
                (*env)->ReleaseStringUTFChars(env, exception, message);
                JavaDeleteRef(exception);

                data[i].molecule = NULL;
            }
            else
            {
                moleculeArray = (jbyteArray) (*env)->GetObjectField(env, resultElement, lucyLoaderMoleculeField);

                jsize moleculeSize = (*env)->GetArrayLength(env, moleculeArray);

                molecule = (*env)->GetByteArrayElements(env, moleculeArray, NULL);
                java_check_exception(__func__);

                data[i].molecule = (bytea *) palloc(VARHDRSZ + moleculeSize);
                SET_VARSIZE(data[i].molecule, VARHDRSZ + moleculeSize);
                memcpy(VARDATA(data[i].molecule), molecule, moleculeSize);

                data[i].error = NULL;

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
        JavaDeleteByteArray(moleculeArray, molecule, JNI_ABORT);

        PG_RE_THROW();
    }
    PG_END_TRY();
}
