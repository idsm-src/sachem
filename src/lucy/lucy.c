#include <postgres.h>
#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <Clownfish/String.h>
#include <Clownfish/Vector.h>
#include <Lucy/Analysis/RegexTokenizer.h>
#include <Lucy/Document/Doc.h>
#include <Lucy/Document/HitDoc.h>
#include <Lucy/Index/Indexer.h>
#include <Lucy/Object/BitVector.h>
#include <Lucy/Plan/FullTextType.h>
#include <Lucy/Plan/Schema.h>
#include <Lucy/Plan/StringType.h>
#include <Lucy/Search/ANDQuery.h>
#include <Lucy/Search/Collector.h>
#include <Lucy/Search/IndexSearcher.h>
#include <Lucy/Search/MatchAllQuery.h>
#include <Lucy/Search/QueryParser.h>
#include <Lucy/Search/TermQuery.h>
#include "lucy.h"


#define ID_NAME     "id"
#define FP_NAME     "fp"
#define FP_PATTERN  "[a-zA-Z0-9+/]+"
#define BOOLOP      "AND"

#define safeDecref(x)           do { void *obj = x; x = NULL; DECREF(obj); } while(0)
#define safeNothrowDecref(x)    do { void *obj = x; x = NULL; nothrowDecref(obj); } while(0)


typedef struct
{
    Lucy *lucy;
    StringType *stype;
    RegexTokenizer *rt;
    String *re;
    FullTextType *fttype;
    String *boolop;
} InitRoutineContext;


typedef struct
{
    Lucy *lucy;
    const char *indexPath;
} SetFolderRoutineContext;


typedef struct
{
    Lucy *lucy;
    int32_t id;
    StringFingerprint fp;
    Doc *doc;
    String *idValue;
    String *fpValue;
} AddRoutineContext;


typedef struct
{
    Lucy *lucy;
    const char *indexPath;
    String *folder;
} AddIndexRoutineContext;


typedef struct
{
    Lucy *lucy;
    int32_t id;
    String *idValue;
} DeleteRoutineContext;


typedef struct
{
    Lucy *lucy;
    StringFingerprint fp;
    int max_results;
    String *queryStr;
    Query *query;
} SearchRoutineContext;


typedef struct
{
    Lucy *lucy;
    LucyResultSet *resultSet;
    int *results;
    size_t size;
    size_t loaded;
    HitDoc *hit;
    String *id;
} GetRoutineContext;


typedef struct
{
    Err *error;
    String *message;
} GetErrorRoutineContext;


static void decref(void *obj)
{
    DECREF(obj);
}


static void nothrowDecref(void *obj)
{
    Err *error = Err_trap((Err_Attempt_t) &decref, obj);

    if(error != NULL)
    {
        safeNothrowDecref(error);
        elog(NOTICE, "cannot destroy a lucy object in error handler");
    }
}


static void throwError(Err *error)
{
    safeNothrowDecref(error);
    elog(ERROR, "lucy error");
}


static void base_init(InitRoutineContext *context)
{
    Lucy *lucy = context->lucy;

    lucy_bootstrap_parcel();

    lucy->idF = Str_new_wrap_trusted_utf8(ID_NAME, sizeof(ID_NAME) - 1);
    lucy->fpF = Str_new_wrap_trusted_utf8(FP_NAME, sizeof(FP_NAME) - 1);

    lucy->schema = Schema_new();
    context->stype = StringType_new();
    StringType_Set_Indexed(context->stype, true);
    Schema_Spec_Field(lucy->schema, lucy->idF, (FieldType *) context->stype);
    safeDecref(context->stype);

    context->re = Str_new_wrap_trusted_utf8(FP_PATTERN, sizeof(FP_PATTERN) - 1);
    context->rt = RegexTokenizer_new(context->re);
    context->fttype = FullTextType_new((Analyzer *) context->rt);
    FullTextType_Set_Indexed(context->fttype, true);
    FullTextType_Set_Highlightable(context->fttype, false);
    FullTextType_Set_Stored(context->fttype, false);
    Schema_Spec_Field(lucy->schema, lucy->fpF, (FieldType *) context->fttype);
    safeDecref(context->re);
    safeDecref(context->rt);
    safeDecref(context->fttype);

    context->boolop = Str_new_wrap_trusted_utf8(BOOLOP, sizeof(BOOLOP) - 1);
    lucy->qparser = QParser_new(lucy->schema, NULL, context->boolop, NULL);
    safeDecref(context->boolop);
}


void lucy_init(Lucy *lucy)
{
    lucy->folder = NULL;
    lucy->idF = NULL;
    lucy->fpF = NULL;
    lucy->schema = NULL;
    lucy->qparser = NULL;
    lucy->indexer = NULL;
    lucy->searcher = NULL;
    lucy->collector = NULL;
    lucy->hits = NULL;

    InitRoutineContext context;
    context.lucy = lucy;
    context.stype = NULL;
    context.rt = NULL;
    context.re = NULL;
    context.fttype = NULL;
    context.boolop = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_init, (void *) &context);

    if(error != NULL)
    {
        safeNothrowDecref(lucy->folder);
        safeNothrowDecref(lucy->idF);
        safeNothrowDecref(lucy->fpF);
        safeNothrowDecref(lucy->schema);
        safeNothrowDecref(lucy->qparser);

        safeNothrowDecref(context.stype);
        safeNothrowDecref(context.rt);
        safeNothrowDecref(context.re);
        safeNothrowDecref(context.fttype);
        safeNothrowDecref(context.boolop);

        throwError(error);
    }
}


static void base_set_folder(SetFolderRoutineContext *context)
{
    Lucy *lucy = context->lucy;

    safeDecref(lucy->folder);
    safeDecref(context->lucy->searcher);
    safeDecref(context->lucy->collector);
    safeDecref(context->lucy->hits);

    lucy->folder = Str_new_from_trusted_utf8(context->indexPath, strlen(context->indexPath));
}


void lucy_set_folder(Lucy *lucy, const char *indexPath)
{
    SetFolderRoutineContext context;
    context.lucy = lucy;
    context.indexPath = indexPath;

    Err *error = Err_trap((Err_Attempt_t) &base_set_folder, (void *) &context);

    if(error != NULL)
    {
        safeNothrowDecref(lucy->folder);
        safeNothrowDecref(lucy->searcher);
        safeNothrowDecref(lucy->collector);
        safeNothrowDecref(lucy->hits);

        throwError(error);
    }
}


static void base_begin(Lucy *lucy)
{
    lucy->indexer = Indexer_new(lucy->schema, (Obj *) (lucy->folder), NULL, Indexer_CREATE);
}


void lucy_begin(Lucy *lucy)
{
    lucy->indexer = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_begin, (void *) lucy);

    if(error != NULL)
        throwError(error);
}


static void base_add(AddRoutineContext *context)
{
    context->doc = Doc_new(NULL, 0);

    char buffer[12];
    int size = sprintf(buffer, "%u", context->id);
    context->idValue = Str_new_wrap_trusted_utf8(buffer, size);
    Doc_Store(context->doc, context->lucy->idF, (Obj *) context->idValue);
    safeDecref(context->idValue);

    context->fpValue = Str_new_wrap_trusted_utf8(context->fp.data, context->fp.size);
    Doc_Store(context->doc, context->lucy->fpF, (Obj *) context->fpValue);
    safeDecref(context->fpValue);

    Indexer_Add_Doc(context->lucy->indexer, context->doc, 1.0);
    safeDecref(context->doc);
}


void lucy_add(Lucy *lucy, int32_t id, StringFingerprint fp)
{
    AddRoutineContext context;
    context.lucy = lucy;
    context.id = id;
    context.fp = fp;
    context.doc = NULL;
    context.idValue = NULL;
    context.fpValue = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_add, (void *) &context);

    if(error != NULL)
    {
        safeNothrowDecref(context.doc);
        safeNothrowDecref(context.idValue);
        safeNothrowDecref(context.fpValue);

        throwError(error);
    }
}


static void base_add_index(AddIndexRoutineContext *context)
{
    context->folder = Str_new_from_trusted_utf8(context->indexPath, strlen(context->indexPath));
    Indexer_Add_Index(context->lucy->indexer, (Obj *) context->folder);
    safeDecref(context->folder);
}


void lucy_add_index(Lucy *lucy, const char *indexPath)
{
    AddIndexRoutineContext context;
    context.lucy = lucy;
    context.indexPath = indexPath;
    context.folder = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_add_index, (void *) &context);

    if(error != NULL)
    {
        safeNothrowDecref(context.folder);

        throwError(error);
    }
}


static void base_delete(DeleteRoutineContext *context)
{
    char buffer[12];
    int size = sprintf(buffer, "%u", context->id);
    context->idValue = Str_new_wrap_trusted_utf8(buffer, size);

    Indexer_Delete_By_Term(context->lucy->indexer, context->lucy->idF, (Obj *) context->idValue);
    safeDecref(context->idValue);
}


void lucy_delete(Lucy *lucy, int32_t id)
{
    DeleteRoutineContext context;
    context.lucy = lucy;
    context.id = id;
    context.idValue = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_delete, (void *) &context);

    if(error != NULL)
    {
        safeNothrowDecref(context.idValue);

        throwError(error);
    }
}


static void base_optimize(Lucy *lucy)
{
    Indexer_Optimize(lucy->indexer);
}


void lucy_optimize(Lucy *lucy)
{
    Err *error = Err_trap((Err_Attempt_t) &base_optimize, (void *) lucy);

    if(error != NULL)
        throwError(error);
}


static void base_commit(Lucy *lucy)
{
    Indexer_Commit(lucy->indexer);

    safeDecref(lucy->indexer);
}


void lucy_commit(Lucy *lucy)
{
    Err *error = Err_trap((Err_Attempt_t) &base_commit, (void *) lucy);

    if(error != NULL)
    {
        safeNothrowDecref(lucy->indexer);
        throwError(error);
    }
}


void lucy_rollback(Lucy *lucy)
{
    safeNothrowDecref(lucy->indexer);
}


static void base_search(SearchRoutineContext *context)
{
    if(context->lucy->searcher == NULL)
    {
        context->lucy->searcher = IxSearcher_new((Obj *) (context->lucy->folder));

        int32_t maxId = IxSearcher_Doc_Max(context->lucy->searcher);
        context->lucy->hits = BitVec_new(maxId + 1);

        context->lucy->collector = (Collector *) BitColl_new(context->lucy->hits);
    }


    if(context->fp.size != 0)
    {
        context->queryStr = Str_new_wrap_trusted_utf8(context->fp.data, context->fp.size);
        context->query = QParser_Parse(context->lucy->qparser, context->queryStr);
    }
    else
    {
        context->query = (Query *) MatchAllQuery_new();
    }

    BitVec_Clear_All(context->lucy->hits);
    IxSearcher_Collect(context->lucy->searcher, context->query, context->lucy->collector);

    safeDecref(context->query);
    safeDecref(context->queryStr);
}


LucyResultSet lucy_search(Lucy *lucy, StringFingerprint fp, int max_results)
{
    SearchRoutineContext context;
    context.lucy = lucy;
    context.fp = fp;
    context.max_results = max_results;
    context.queryStr = NULL;
    context.query = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_search, (void *) &context);

    if(error != NULL)
    {
        if(lucy->searcher == NULL || lucy->collector == NULL || lucy->hits == NULL)
        {
            safeNothrowDecref(lucy->searcher);
            safeNothrowDecref(lucy->collector);
            safeNothrowDecref(lucy->hits);
        }

        safeNothrowDecref(context.query);
        safeNothrowDecref(context.queryStr);

        throwError(error);
    }

    return (LucyResultSet) { .possition = 0 };
}


static void base_get(GetRoutineContext *context)
{
    int ret = 0;
    int size = context->size;
    int *results = context->results;

    while(size > 0)
    {
        size_t docId = BitVec_Next_Hit(context->lucy->hits, context->resultSet->possition);
        context->resultSet->possition = docId == -1 ? -1 : docId + 1;

        if(docId == -1)
            break;

        context->hit = IxSearcher_Fetch_Doc(context->lucy->searcher, docId);

        context->id = (String *) HitDoc_Extract(context->hit, context->lucy->idF);
        *(results++) = Str_To_I64(context->id);

        ret++;
        size--;

        safeDecref(context->id);
        safeDecref(context->hit);
    }

    context->loaded = ret;
}


size_t lucy_get(Lucy *lucy, LucyResultSet *resultSet, int *results, size_t size)
{
    GetRoutineContext context;
    context.lucy = lucy;
    context.resultSet = resultSet;
    context.results = results;
    context.size = size;
    context.loaded = 0;
    context.hit = NULL;
    context.id = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_get, (void *) &context);

    if(error != NULL)
    {
        safeNothrowDecref(context.hit);
        safeNothrowDecref(context.id);

        throwError(error);
    }

    return context.loaded;
}


void lucy_fail(Lucy *lucy, LucyResultSet *resultSet)
{
}


static void link_directory_at(int olddirfd, int newdirfd)
{
    DIR *dp = NULL;
    int subfd = -1;
    int newsubfd = -1;

    PG_TRY();
    {
        if((dp = fdopendir(olddirfd)) == NULL)
            elog(ERROR, "%s: fdopendir() failed", __func__);


        struct dirent *ep;

        while((ep = readdir(dp)))
        {
            if(ep->d_name[0] == '.')
                continue;

            if((subfd = openat(olddirfd, ep->d_name, 0)) == -1)
                elog(ERROR, "%s: openat() failed", __func__);

            struct stat statbuf;

            if(fstat(subfd, &statbuf) != 0)
                elog(ERROR, "%s: fstat() failed", __func__);

            if(statbuf.st_mode & S_IFDIR)
            {
                if(mkdirat(newdirfd, ep->d_name, 0700) != 0)
                    elog(ERROR, "%s: mkdirat() failed", __func__);

                if((newsubfd = openat(newdirfd, ep->d_name, 0)) == -1)
                    elog(ERROR, "%s: openat() failed", __func__);

                link_directory_at(subfd, newsubfd);

                close(newsubfd);
                newsubfd = -1;
            }
            else if(statbuf.st_mode & S_IFREG)
            {
                if(linkat(olddirfd, ep->d_name, newdirfd, ep->d_name, 0) != 0)
                    elog(ERROR, "%s: linkat() failed", __func__);
            }

            close(subfd);
            subfd = -1;
        }

        closedir(dp);
        dp = NULL;
    }
    PG_CATCH();
    {
        if(newsubfd != -1)
            close(newsubfd);

        if(subfd != -1)
            close(subfd);

        if(dp != NULL)
            closedir(dp);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


void lucy_link_directory(const char *oldPath, const char *newPath)
{
    if(mkdir(newPath, 0700) != 0)
        elog(ERROR, "%s: mkdir() failed", __func__);


    if(oldPath != NULL)
    {
        int newfd = -1;
        int oldfd = -1;

        PG_TRY();
        {
            if((newfd = open(newPath, O_DIRECTORY)) == -1)
                elog(ERROR, "%s: open() failed", __func__);

            if((oldfd = open(oldPath, O_DIRECTORY)) == -1)
                elog(ERROR, "%s: open() failed", __func__);

            link_directory_at(oldfd, newfd);

            close(oldfd);
            oldfd = -1;

            close(newfd);
            newfd = -1;
        }
        PG_CATCH();
        {
            if(oldfd != -1)
                close(oldfd);

            if(newfd != -1)
                close(newfd);

            PG_RE_THROW();
        }
        PG_END_TRY();
    }
}


static void unlink_directory_at(int dirfd)
{
    DIR *dp = NULL;
    int subfd = -1;

    PG_TRY();
    {
        if((dp = fdopendir(dirfd)) == NULL)
            elog(ERROR, "%s: fdopendir() failed", __func__);


        struct dirent *ep;

        while((ep = readdir(dp)))
        {
            if(ep->d_name[0] == '.')
                continue;

            if((subfd = openat(dirfd, ep->d_name, 0)) == -1)
                elog(ERROR, "%s: openat() failed", __func__);

            struct stat statbuf;
            if(fstat(subfd, &statbuf) != 0)
                elog(ERROR, "%s: fstat() failed", __func__);

            if(statbuf.st_mode & S_IFDIR)
            {
                unlink_directory_at(subfd);

                if(unlinkat(dirfd, ep->d_name, AT_REMOVEDIR) != 0)
                   elog(ERROR, "%s: unlinkat() failed", __func__);
            }
            else if(statbuf.st_mode & S_IFREG)
            {
                if(unlinkat(dirfd, ep->d_name, 0) != 0)
                    elog(ERROR, "%s: unlinkat() failed", __func__);
            }

            close(subfd);
            subfd = -1;
        }

        closedir(dp);
        dp = NULL;
    }
    PG_CATCH();
    {
        if(subfd != -1)
            close(subfd);

        if(dp != NULL)
            closedir(dp);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


void lucy_delete_directory(const char *path)
{
    int fd = -1;

    PG_TRY();
    {
        if((fd = open(path, O_DIRECTORY)) == -1)
            elog(ERROR, "%s: open() failed", __func__);

        unlink_directory_at(fd);

        close(fd);
        fd = -1;

        if(rmdir(path) != 0)
           elog(ERROR, "%s: rmdir() failed", __func__);
    }
    PG_CATCH();
    {
        if(fd != -1)
            close(fd);

        PG_RE_THROW();
    }
    PG_END_TRY();
}
