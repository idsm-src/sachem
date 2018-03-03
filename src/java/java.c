#include <postgres.h>
#include <libpq/pqsignal.h>
#include <tcop/tcopprot.h>
#include <pthread.h>
#include "java.h"
#include "sachem.h"


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
            .nOptions = 2,
            .options = (JavaVMOption []) {
                (JavaVMOption) { .optionString = "-Djava.class.path=" JARDIR "/sachem.jar:" JARDIR "/cdk-2.1.1.jar:" JARDIR "/lucene-core-7.2.1.jar"},
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
            elog(ERROR, "%s: cannot initialize JVM", __func__);
    }

    exceptionClass = (*env)->FindClass(env, "java/lang/Throwable");
    java_check_exception(__func__);

    toStringMethod = (*env)->GetMethodID(env, exceptionClass, "toString", "()Ljava/lang/String;");
    java_check_exception(__func__);


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

        const char *mstr = (*env)->GetStringUTFChars(env, message, NULL);
        error = (char *) palloc(strlen(mstr) + 1);
        strcpy(error, mstr);

        (*env)->ReleaseStringUTFChars(env, message, mstr);
        JavaDeleteRef(message);
        JavaDeleteRef(exception);
    }

    elog(ERROR, "%s: java error: %s", str, error);
}


void __attribute__ ((destructor)) java_finish(void)
{
    if(jvm != NULL)
        (*jvm)->DestroyJavaVM(jvm);
}
