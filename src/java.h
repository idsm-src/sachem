#ifndef JAVA_H_
#define JAVA_H_

#include <postgres.h>
#include <utils/array.h>
#include <jni.h>
#include <stdbool.h>
#include <stdint.h>


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


#define JavaDeleteIntegerArray(array, buffer, mode) \
    do { \
        if(likely(array != NULL && buffer != NULL)) (*env)->ReleaseIntArrayElements(env, array, buffer, mode); \
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



typedef struct
{
    uint8_t *molecule;
    bool *restH;
} SubstructureQueryData;


typedef struct
{
    int16_t *counts;
    int16_t *fp;
    uint8_t *molecule;
    bool *restH;

    int fpLength;
} OrchemSubstructureQueryData;


typedef struct
{
    int bitCount;
    ArrayType *fp;
    ArrayType *counts;
    bytea *molecule;
    text *error;
} OrchemLoaderData;


typedef struct
{
    bytea *molecule;
    text *error;
} LucyLoaderData;


extern JavaVM* jvm;
extern JNIEnv* env;


void java_check_exception(const char *str);
int java_parse_substructure_query(SubstructureQueryData **data, char* query, size_t queryLength, int32_t type, bool implicitHydrogens, bool tautomers);
int java_parse_orchem_substructure_query(OrchemSubstructureQueryData **data, char* query, size_t queryLength, int32_t type, bool implicitHydrogens, bool tautomers);
int java_parse_orchem_similarity_query(uint64_t **data, char* query, size_t queryLength, int32_t type);
void java_parse_orchem_data(size_t count, VarChar **molfiles, OrchemLoaderData *data);
void java_parse_lucy_data(size_t count, VarChar **molfiles, LucyLoaderData *data);
void java_module_init(void);
void java_module_finish(void);

#endif /* JAVA_H_ */
