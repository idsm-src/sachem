#include <postgres.h>
#include <stdio.h>
#include <Clownfish/String.h>
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
#include <Lucy/Search/QueryParser.h>
#include <Lucy/Search/TermQuery.h>
#include "lucy.h"


#define ID_NAME     "id"
#define FP_NAME     "fp"
#define FP_PATTERN  "[a-zA-Z0-9+/]+"
#define BOOLOP      "AND"

#define saveDecref(x)           do { void *obj = x; x = NULL; DECREF(obj); } while(0)
#define saveNothrowDecref(x)    do { void *obj = x; x = NULL; nothrowDecref(obj); } while(0)


typedef struct
{
    Lucy *lucy;
    const char *indexPath;
    StringType *stype;
    RegexTokenizer *rt;
    String *re;
    FullTextType *fttype;
    String *boolop;
} InitRoutineContext;


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
        saveNothrowDecref(error);
        elog(NOTICE, "cannot destroy a lucy object in error handler");
    }
}


static void throwError(Err *error)
{
    saveNothrowDecref(error);
    elog(ERROR, "lucy error");
}


static void base_init(InitRoutineContext *context)
{
    Lucy *lucy = context->lucy;

    lucy_bootstrap_parcel();

    lucy->folder = Str_new_from_trusted_utf8(context->indexPath, strlen(context->indexPath));
    lucy->idF = Str_new_wrap_trusted_utf8(ID_NAME, sizeof(ID_NAME) - 1);
    lucy->fpF = Str_new_wrap_trusted_utf8(FP_NAME, sizeof(FP_NAME) - 1);

    lucy->schema = Schema_new();
    context->stype = StringType_new();
    StringType_Set_Indexed(context->stype, false);
    Schema_Spec_Field(lucy->schema, lucy->idF, (FieldType *) context->stype);
    saveDecref(context->stype);

    context->re = Str_new_wrap_trusted_utf8(FP_PATTERN, sizeof(FP_PATTERN) - 1);
    context->rt = RegexTokenizer_new(context->re);
    context->fttype = FullTextType_new((Analyzer *) context->rt);
    FullTextType_Set_Indexed(context->fttype, true);
    FullTextType_Set_Highlightable(context->fttype, false);
    FullTextType_Set_Stored(context->fttype, false);
    Schema_Spec_Field(lucy->schema, lucy->fpF, (FieldType *) context->fttype);
    saveDecref(context->re);
    saveDecref(context->rt);
    saveDecref(context->fttype);

    context->boolop = Str_new_wrap_trusted_utf8(BOOLOP, sizeof(BOOLOP) - 1);
    lucy->qparser = QParser_new(lucy->schema, NULL, context->boolop, NULL);
    saveDecref(context->boolop);
}


void lucy_init(Lucy *lucy, const char *indexPath)
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
    context.indexPath = indexPath;
    context.stype = NULL;
    context.rt = NULL;
    context.re = NULL;
    context.fttype = NULL;
    context.boolop = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_init, (void *) &context);

    if(error != NULL)
    {
        saveNothrowDecref(lucy->folder);
        saveNothrowDecref(lucy->idF);
        saveNothrowDecref(lucy->fpF);
        saveNothrowDecref(lucy->schema);
        saveNothrowDecref(lucy->qparser);

        saveNothrowDecref(context.stype);
        saveNothrowDecref(context.rt);
        saveNothrowDecref(context.re);
        saveNothrowDecref(context.fttype);
        saveNothrowDecref(context.boolop);

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
    saveDecref(context->idValue);

    context->fpValue = Str_new_wrap_trusted_utf8(context->fp.data, context->fp.size);
    Doc_Store(context->doc, context->lucy->fpF, (Obj *) context->fpValue);
    saveDecref(context->fpValue);

    Indexer_Add_Doc(context->lucy->indexer, context->doc, 1.0);
    saveDecref(context->doc);
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
        saveNothrowDecref(context.doc);
        saveNothrowDecref(context.idValue);
        saveNothrowDecref(context.fpValue);

        throwError(error);
    }
}


static void base_delete(DeleteRoutineContext *context)
{
    char buffer[12];
    int size = sprintf(buffer, "%u", context->id);
    context->idValue = Str_new_wrap_trusted_utf8(buffer, size);

    Indexer_Delete_By_Term(context->lucy->indexer, context->lucy->idF, (Obj *) context->idValue);
    saveDecref(context->idValue);
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
        saveNothrowDecref(context.idValue);

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
    saveDecref(lucy->indexer);
}


void lucy_commit(Lucy *lucy)
{
    Err *error = Err_trap((Err_Attempt_t) &base_commit, (void *) lucy);

    if(error != NULL)
        throwError(error);
}


void lucy_rollback(Lucy *lucy)
{
    saveNothrowDecref(lucy->indexer);
}


static void base_search(SearchRoutineContext *context)
{
    if(context->lucy->searcher == NULL)
        context->lucy->searcher = IxSearcher_new((Obj *) (context->lucy->folder));

    context->queryStr = Str_new_wrap_trusted_utf8(context->fp.data, context->fp.size);
    context->query = QParser_Parse(context->lucy->qparser, context->queryStr);
    context->hits = IxSearcher_Hits(context->lucy->searcher, (Obj *) context->query, 0, context->max_results, NULL);

    saveDecref(context->query);
    saveDecref(context->queryStr);
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
        saveNothrowDecref(context.lucy->searcher);
        saveNothrowDecref(context.query);
        saveNothrowDecref(context.queryStr);
        saveNothrowDecref(context.hits);

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
            saveDecref(context->resultSet->hits);
            break;
        }

        context->id = (String *) HitDoc_Extract(context->hit, context->lucy->idF);
        *(results++) = Str_To_I64(context->id);

        ret++;
        size--;

        saveDecref(context->id);
        saveDecref(context->hit);
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
        saveNothrowDecref(context.hit);
        saveNothrowDecref(context.id);

        throwError(error);
    }

    return context.loaded;
}


void lucy_fail(Lucy *lucy, LucyResultSet *resultSet)
{
    saveNothrowDecref(resultSet->hits);
}
