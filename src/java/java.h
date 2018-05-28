#ifndef JAVA_H_
#define JAVA_H_

#include <jni.h>


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


extern JavaVM* jvm;
extern JNIEnv* env;


void java_check_exception(const char *str);
void java_init(void);
void java_terminate(void);

#endif /* JAVA_H_ */
