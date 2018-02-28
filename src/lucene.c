#include <postgres.h>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "sachem.h"
#include "lucene.h"
#include "java.h"


static jclass luceneClass;
static jmethodID setFolderMethod;
static jmethodID beginMethod;
static jmethodID addMethod;
static jmethodID deleteMethod;
static jmethodID optimizeMethod;
static jmethodID commitMethod;
static jmethodID rollbackMethod;
static jmethodID searchMethod;
static jmethodID getMethod;


void lucene_init(Lucene *lucene)
{
    luceneClass = (*env)->FindClass(env, "cz/iocb/sachem/lucene/Lucene");

    setFolderMethod = (*env)->GetMethodID(env, luceneClass, "setFolder", "(Ljava/lang/String;)V");
    java_check_exception(__func__);

    beginMethod = (*env)->GetMethodID(env, luceneClass, "begin", "()V");
    java_check_exception(__func__);

    addMethod = (*env)->GetMethodID(env, luceneClass, "add", "(I[I)V");
    java_check_exception(__func__);

    deleteMethod = (*env)->GetMethodID(env, luceneClass, "delete", "(I)V");
    java_check_exception(__func__);

    optimizeMethod = (*env)->GetMethodID(env, luceneClass, "optimize", "()V");
    java_check_exception(__func__);

    commitMethod = (*env)->GetMethodID(env, luceneClass, "commit", "()V");
    java_check_exception(__func__);

    rollbackMethod = (*env)->GetMethodID(env, luceneClass, "rollback", "()V");
    java_check_exception(__func__);

    searchMethod = (*env)->GetMethodID(env, luceneClass, "search", "([II)Lcz/iocb/sachem/lucene/Lucene$ResultSet;");
    java_check_exception(__func__);

    getMethod = (*env)->GetMethodID(env, luceneClass, "get", "(Lcz/iocb/sachem/lucene/Lucene$ResultSet;[I)I");
    java_check_exception(__func__);

    jmethodID constructor = (*env)->GetMethodID(env, luceneClass, "<init>", "()V");
    java_check_exception(__func__);

    *lucene = (*env)->NewObject(env, luceneClass, constructor);
    java_check_exception(__func__);
}


void lucene_set_folder(Lucene *lucene, const char *path)
{
    jstring folder = NULL;

    PG_TRY();
    {
        folder = (*env)->NewStringUTF(env, path);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, *lucene, setFolderMethod, folder);
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


void lucene_begin(Lucene *lucene)
{
    (*env)->CallVoidMethod(env, *lucene, beginMethod);
    java_check_exception(__func__);
}


void lucene_add(Lucene *lucene, int32_t id, IntegerFingerprint fp)
{
    jintArray fpArray = NULL;

    PG_TRY();
    {
        fpArray = (jintArray) (*env)->NewIntArray(env, fp.size);
        java_check_exception(__func__);

        (*env)->SetIntArrayRegion(env, fpArray, 0, fp.size, (jint*) fp.data);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, *lucene, addMethod, (jint) id, fpArray);
        java_check_exception(__func__);

        JavaDeleteRef(fpArray);
    }
    PG_CATCH();
    {
        JavaDeleteRef(fpArray);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


void lucene_delete(Lucene *lucene, int32_t id)
{
    (*env)->CallVoidMethod(env, *lucene, deleteMethod, (jint) id);
    java_check_exception(__func__);
}


void lucene_optimize(Lucene *lucene)
{
    (*env)->CallVoidMethod(env, *lucene, optimizeMethod);
    java_check_exception(__func__);
}


void lucene_commit(Lucene *lucene)
{
    (*env)->CallVoidMethod(env, *lucene, commitMethod);
    java_check_exception(__func__);
}


void lucene_rollback(Lucene *lucene)
{
    (*env)->CallVoidMethod(env, *lucene, rollbackMethod);
    java_check_exception(__func__);
}


LuceneResultSet lucene_search(Lucene *lucene, IntegerFingerprint fp, int maxResultCount)
{
    jobject result = NULL;
    jintArray fpArray = NULL;

    PG_TRY();
    {
        fpArray = (jintArray) (*env)->NewIntArray(env, fp.size);
        java_check_exception(__func__);

        (*env)->SetIntArrayRegion(env, fpArray, 0, fp.size, (jint*) fp.data);
        java_check_exception(__func__);

        result = (*env)->CallObjectMethod(env, *lucene, searchMethod, fpArray, (jint) maxResultCount);
        java_check_exception(__func__);

        JavaDeleteRef(fpArray);
    }
    PG_CATCH();
    {
        JavaDeleteRef(fpArray);
        JavaDeleteRef(result);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return result;
}


size_t lucene_get(Lucene *lucene, LuceneResultSet *resultSet, int *buffer, size_t size)
{
    jintArray resultArray = NULL;
    jint *results = NULL;
    size_t count = 0;

    PG_TRY();
    {
        resultArray = (jintArray) (*env)->NewIntArray(env, size);
        java_check_exception(__func__);

        count = (*env)->CallIntMethod(env, *lucene, getMethod, *resultSet, resultArray);
        java_check_exception(__func__);

        results = (*env)->GetIntArrayElements(env, resultArray, 0);
        java_check_exception(__func__);

        for(int i = 0; i < count; i++)
            buffer[i] = results[i];

        if(count == 0)
            JavaDeleteRef(*resultSet);

        JavaDeleteIntegerArray(resultArray, results, JNI_ABORT);
    }
    PG_CATCH();
    {
        JavaDeleteIntegerArray(resultArray, results, JNI_ABORT);
        JavaDeleteRef(*resultSet);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return count;
}


void lucene_fail(Lucene *lucene, LuceneResultSet *resultSet)
{
    JavaDeleteRef(*resultSet);
}


void lucene_link_directory(const char *oldPath, const char *newPath)
{
    int newdirfd = -1;
    int olddirfd = -1;
    DIR *dp = NULL;


    PG_TRY();
    {
        if(mkdir(newPath, 0700) != 0)
            elog(ERROR, "%s: mkdir() failed", __func__);


        if(oldPath != NULL)
        {
            if((newdirfd = open(newPath, 0)) == -1)
                elog(ERROR, "%s: open() failed", __func__);

            if((olddirfd = open(oldPath, 0)) == -1)
                elog(ERROR, "%s: open() failed", __func__);

            if((dp = fdopendir(olddirfd)) == NULL)
                elog(ERROR, "%s: fdopendir() failed", __func__);


            struct dirent *ep;

            while(ep = readdir (dp))
            {
                if(ep->d_name[0] == '.')
                    continue;

                if(linkat(olddirfd, ep->d_name, newdirfd, ep->d_name, 0) != 0)
                    elog(ERROR, "%s: linkat() failed", __func__);
            }

            closedir(dp);
            dp = NULL;

            close(olddirfd);
            olddirfd = -1;

            close(newdirfd);
            newdirfd = -1;
        }
    }
    PG_CATCH();
    {
        if(dp != NULL)
            closedir(dp);

        if(olddirfd != -1)
            close(olddirfd);

        if(newdirfd != -1)
            close(newdirfd);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


void lucene_delete_directory(const char *path)
{
    int dirfd = -1;
    DIR *dp = NULL;

    PG_TRY();
    {
        if((dirfd = open(path, 0)) == -1)
            elog(ERROR, "%s: open() failed", __func__);

        if((dp = fdopendir(dirfd)) == NULL)
            elog(ERROR, "%s: fdopendir() failed", __func__);


        struct dirent *ep;

        while(ep = readdir(dp))
        {
            if(ep->d_name[0] == '.')
                continue;

            if(unlinkat(dirfd, ep->d_name, 0) != 0)
                elog(ERROR, "%s: unlinkat() failed", __func__);
        }

        closedir(dp);
        dp = NULL;

        close(dirfd);
        dirfd = -1;

        if(rmdir(path) != 0)
           elog(ERROR, "%s: rmdir() failed", __func__);
    }
    PG_CATCH();
    {
        if(dp != NULL)
            closedir(dp);

        if(dirfd != -1)
            close(dirfd);

        PG_RE_THROW();
    }
    PG_END_TRY();
}
