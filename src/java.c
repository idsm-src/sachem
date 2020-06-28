#include <postgres.h>
#include <pthread.h>
#include <libpq/pqsignal.h>
#include <tcop/tcopprot.h>
#include "java.h"
#include "inchi.h"
#include "isomorphism.h"


JavaVM* jvm = NULL;
JNIEnv* env = NULL;

static bool initialised = false;
static pthread_t mainThread;
static jclass exceptionClass = NULL;
static jmethodID toStringMethod = NULL;


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


void java_init(void)
{
    if(likely(initialised))
        return;

    if(jvm == NULL || env == NULL)
    {
        JavaVMInitArgs args = (JavaVMInitArgs) {
            .version = JNI_VERSION_1_8,
            .nOptions = 1,
            .options = (JavaVMOption []) { (JavaVMOption) {
                .optionString = "-Djava.class.path=" JARDIR "/sachem.jar:"
                                                     JARDIR "/beam-core-1.3.3.jar:"
                                                     JARDIR "/cdk-atomtype-2.3.jar:"
                                                     JARDIR "/cdk-core-2.3.jar:"
                                                     JARDIR "/cdk-ctab-2.3.jar:"
                                                     JARDIR "/cdk-interfaces-2.3.jar:"
                                                     JARDIR "/cdk-ioformats-2.3.jar:"
                                                     JARDIR "/cdk-isomorphism-2.3.jar:"
                                                     JARDIR "/cdk-silent-2.3.jar:"
                                                     JARDIR "/cdk-smiles-2.3.jar:"
                                                     JARDIR "/cdk-standard-2.3.jar:"
                                                     JARDIR "/cdk-valencycheck-2.3.jar:"
                                                     JARDIR "/guava-29.0-jre.jar:"
                                                     JARDIR "/log4j-1.2-api-2.13.3.jar:"
                                                     JARDIR "/log4j-api-2.13.3.jar:"
                                                     JARDIR "/log4j-core-2.13.3.jar:"
                                                     JARDIR "/lucene-core-8.5.2.jar:"
                                                     JARDIR "/vecmath-1.5.2.jar" }},
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
    }

    exceptionClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "java/lang/Throwable"));
    java_check_exception(__func__);

    toStringMethod = (*env)->GetMethodID(env, exceptionClass, "toString", "()Ljava/lang/String;");
    java_check_exception(__func__);

    isomorphism_init();
    inchi_init();

    initialised = true;
}


void java_check_exception(const char *str)
{
    jthrowable exception = (*env)->ExceptionOccurred(env);

    if(exception == NULL)
        return;

    (*env)->ExceptionDescribe(env);
    (*env)->ExceptionClear(env);
    char *error = NULL;

    if(toStringMethod != NULL)
    {
        jstring message = (jstring)(*env)->CallObjectMethod(env, exception, toStringMethod);

        const char *mstr = message != NULL ? (*env)->GetStringUTFChars(env, message, NULL) : NULL;

        if(mstr != NULL)
        {
            error = (char *) palloc(strlen(mstr) + 1);
            strcpy(error, mstr);
            (*env)->ReleaseStringUTFChars(env, message, mstr);
        }
        else
        {
            error = "unknown jvm error";
        }

        JavaDeleteRef(message);
    }

    JavaDeleteRef(exception);

    elog(ERROR, "%s: java error: %s", str, error);
}


static void __attribute__ ((destructor)) java_finish(void)
{
    if(jvm != NULL)
        (*jvm)->DestroyJavaVM(jvm);
}
