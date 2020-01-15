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


extern JNIEnv* env;


void java_check_exception(const char *str);
void java_init(void);

#endif /* JAVA_H_ */
