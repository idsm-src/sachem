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
#include <Lucy/Plan/FullTextType.h>
#include <Lucy/Plan/Schema.h>
#include <Lucy/Plan/StringType.h>
#include <Lucy/Search/ANDQuery.h>
#include <Lucy/Search/Hits.h>
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
    Fingerprint fp;
    Doc *doc;
    String *idValue;
    String *fpValue;
} AddRoutineContext;


typedef struct
{
    Lucy *lucy;
    int32_t id;
    String *idValue;
} DeleteRoutineContext;


typedef struct
{
    Lucy *lucy;
    Fingerprint fp;
    int max_results;
    String *queryStr;
    Query *query;
    Hits *hits;
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
    StringType_Set_Indexed(context->stype, false);
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

    if(lucy->folder != NULL)
        safeDecref(lucy->folder);

    if(context->lucy->searcher != NULL)
        safeDecref(context->lucy->searcher);

    lucy->folder = Str_new_from_trusted_utf8(context->indexPath, strlen(context->indexPath));
}


void lucy_set_folder(Lucy *lucy, const char *indexPath)
{
    SetFolderRoutineContext context;
    context.lucy = lucy;
    context.indexPath = indexPath;

    Err *error = Err_trap((Err_Attempt_t) &base_set_folder, (void *) &context);

    if(error != NULL)
        throwError(error);
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


void lucy_add(Lucy *lucy, int32_t id, Fingerprint fp)
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
        context->lucy->searcher = IxSearcher_new((Obj *) (context->lucy->folder));


    if(context->fp.size != 0)
    {
        context->queryStr = Str_new_wrap_trusted_utf8(context->fp.data, context->fp.size);
        context->query = QParser_Parse(context->lucy->qparser, context->queryStr);
    }
    else
    {
        context->query = (Query *) MatchAllQuery_new();
    }

    context->hits = IxSearcher_Hits(context->lucy->searcher, (Obj *) context->query, 0, context->max_results, NULL);

    safeDecref(context->query);
    safeDecref(context->queryStr);
}


LucyResultSet lucy_search(Lucy *lucy, Fingerprint fp, int max_results)
{
    SearchRoutineContext context;
    context.lucy = lucy;
    context.fp = fp;
    context.max_results = max_results;
    context.queryStr = NULL;
    context.query = NULL;
    context.hits = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_search, (void *) &context);

    if(error != NULL)
    {
        safeNothrowDecref(context.query);
        safeNothrowDecref(context.queryStr);
        safeNothrowDecref(context.hits);

        throwError(error);
    }

    return (LucyResultSet) { .hits = context.hits };
}


static void base_get(GetRoutineContext *context)
{
    int ret = 0;
    int size = context->size;
    int *results = context->results;

    while(size > 0)
    {
        context->hit = Hits_Next(context->resultSet->hits);

        if(context->hit == NULL)
        {
            safeDecref(context->resultSet->hits);
            break;
        }

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
    safeNothrowDecref(resultSet->hits);
}


static void link_directory_at(int olddirfd, int newdirfd)
{
    DIR *dp = fdopendir(olddirfd);

    if(dp == NULL)
        elog(ERROR, "%s: fdopendir() failed", __func__);


    struct dirent *ep;

    while(ep = readdir (dp))
    {
        if(ep->d_name[0] == '.')
            continue;

        int subfd = openat(olddirfd, ep->d_name, 0);

        if(subfd == -1)
            elog(ERROR, "%s: openat() failed", __func__);

        struct stat statbuf;
        if(fstat(subfd, &statbuf) != 0)
            elog(ERROR, "%s: fstat() failed", __func__);

        if(statbuf.st_mode & S_IFDIR)
        {
            if(mkdirat(newdirfd, ep->d_name, 0700) != 0)
                elog(ERROR, "%s: mkdirat() failed", __func__);

            int newsubfd = openat(newdirfd, ep->d_name, 0);

            if(newsubfd == -1)
                elog(ERROR, "%s: openat() failed", __func__);

            link_directory_at(subfd, newsubfd);
        }
        else if(statbuf.st_mode & S_IFREG)
        {
             if(linkat(olddirfd, ep->d_name, newdirfd, ep->d_name, 0) != 0)
                elog(ERROR, "%s: linkat() failed", __func__);
        }
    }

    (void) closedir (dp);
}


void lucy_link_directory(const char *oldPath, const char *newPath)
{
    if(mkdir(newPath, 0700) != 0)
        elog(ERROR, "%s: mkdir() failed", __func__);


    if(oldPath != NULL)
    {
        int newfd = open(newPath, 0);

        if(newfd == -1)
            elog(ERROR, "%s: open() failed", __func__);

        int oldfd = open(oldPath, 0);

        if(oldfd == -1)
            elog(ERROR, "%s: open() failed", __func__);

        link_directory_at(oldfd, newfd);
    }
}
