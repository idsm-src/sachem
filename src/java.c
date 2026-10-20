#include <postgres.h>
#include <pthread.h>
#include <signal.h>
#include <miscadmin.h>
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
            .version = JNI_VERSION_10,
            .nOptions = 3,
            .options = (JavaVMOption []) {
                (JavaVMOption) { .optionString = "--enable-native-access=ALL-UNNAMED" },
                (JavaVMOption) { .optionString = "--add-modules=jdk.incubator.vector" },
                (JavaVMOption) { .optionString = "-Djava.class.path="
                        JARDIR "/sachem.jar:"
                        JARDIR "/beam-core-1.3.10.jar:"
                        JARDIR "/cdk-atomtype-2.13.jar:"
                        JARDIR "/cdk-core-2.13.jar:"
                        JARDIR "/cdk-ctab-2.13.jar:"
                        JARDIR "/cdk-interfaces-2.13.jar:"
                        JARDIR "/cdk-ioformats-2.13.jar:"
                        JARDIR "/cdk-isomorphism-2.13.jar:"
                        JARDIR "/cdk-silent-2.13.jar:"
                        JARDIR "/cdk-smiles-2.13.jar:"
                        JARDIR "/cdk-standard-2.13.jar:"
                        JARDIR "/cdk-valencycheck-2.13.jar:"
                        JARDIR "/guava-33.7.1-jre.jar:"
                        JARDIR "/log4j-1.2-api-2.26.1.jar:"
                        JARDIR "/log4j-api-2.26.1.jar:"
                        JARDIR "/log4j-core-2.26.1.jar:"
                        JARDIR "/lucene-core-10.5.1.jar:"
                        JARDIR "/vecmath-1.5.2.jar" }},
            .ignoreUnrecognized = JNI_FALSE
        };

        /*
         * The signals are ignored while the JVM starts, so that it does not install its own handlers for them, and
         * blocked, so that a cancel or termination request arriving meanwhile is not lost: Linux keeps a blocked signal
         * pending even if it is ignored. Ignoring discards a signal that is already pending, so such a signal is raised
         * again afterwards.
         */
        sigset_t signals;
        sigset_t pending;
        sigset_t mask;

        sigemptyset(&signals);
        sigaddset(&signals, SIGINT);
        sigaddset(&signals, SIGTERM);
        sigaddset(&signals, SIGQUIT);

        pthread_sigmask(SIG_BLOCK, &signals, &mask);
        sigpending(&pending);

        pqsignal(SIGINT, SIG_IGN);
        pqsignal(SIGTERM, SIG_IGN);
        pqsignal(SIGQUIT, SIG_IGN);

        JNI_CreateJavaVM(&(jvm), (void **) &(env), &args);

        mainThread = pthread_self();
        pqsignal(SIGINT,  java_statement_cancel_handler);
        pqsignal(SIGTERM, java_die_handler);
        pqsignal(SIGQUIT, java_quick_die_handler);

        for(int signum = 1; signum < NSIG; signum++)
            if(sigismember(&signals, signum) == 1 && sigismember(&pending, signum) == 1)
                raise(signum);

        pthread_sigmask(SIG_SETMASK, &mask, NULL);

        if(jvm == NULL || env == NULL)
            elog(ERROR, "cannot initialize JVM");

        /* the query is cancelled now if it was cancelled while the JVM started */
        CHECK_FOR_INTERRUPTS();
    }

    exceptionClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "java/lang/Throwable"));
    java_check_exception(__func__);

    toStringMethod = (*env)->GetMethodID(env, exceptionClass, "toString", "()Ljava/lang/String;");
    java_check_exception(__func__);

    native_isomorphism_init(env);
    java_check_exception(__func__);

    inchi_init(env);
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
