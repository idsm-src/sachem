#include <postgres.h>
#include <stdio.h>
#include <Clownfish/String.h>
#include <Clownfish/Vector.h>
#include <Lucy/Analysis/RegexTokenizer.h>
#include <Lucy/Document/Doc.h>
#include <Lucy/Document/HitDoc.h>
#include <Lucy/Index/Indexer.h>
#include <Lucy/Index/Snapshot.h>
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
    bytea *binarySnapshot;
} CommitRoutineContext;


typedef struct
{
    Lucy *lucy;
    bytea *binarySnapshot;
    Snapshot *snapshot;
}
SetSnapshotRoutineContext;


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


static bytea *snapshot_export(Snapshot *snapshot)
{
    Vector *list = Snapshot_List(snapshot);
    uint32_t count = Snapshot_Num_Entries(snapshot);
    size_t size = VARHDRSZ;

    String *path = Snapshot_Get_Path(snapshot);
    size_t pathLength = Str_Get_Size(path);
    size += sizeof(size_t) + pathLength;

    for(int i = 0; i < count; i++)
    {
        String *entry = (String *) Vec_Fetch(list, i);
        size += sizeof(size_t) + Str_Get_Size(entry);
    }

    bytea *result = palloc_extended(size, MCXT_ALLOC_NO_OOM);

    if(result == NULL)
    {
        DECREF(list);
        CFISH_THROW(CFISH_ERR, "not enough memory");
    }

    SET_VARSIZE(result, size);

    char *data = ((char *) result) + VARHDRSZ;

    memcpy(data, &pathLength, sizeof(size_t));
    data += sizeof(size_t);

    memcpy(data, Str_Get_Ptr8(path), pathLength);
    data += pathLength;

    for(int i = 0; i < count; i++)
    {
        String *entry = (String *) Vec_Fetch(list, i);
        size_t length = Str_Get_Size(entry);

        memcpy(data, &length, sizeof(size_t));
        data += sizeof(size_t);

        memcpy(data, Str_Get_Ptr8(entry), length);
        data += length;
    }

    DECREF(list);
    return result;
}


static Snapshot *snapshot_import(bytea *data)
{
    int size = VARSIZE(data) - VARHDRSZ;

    char *buffer = VARDATA(data);
    char *end = buffer + size;
    Snapshot *snapshot = Snapshot_new();

    size_t length = *((size_t *) buffer);
    buffer += sizeof(size_t);

    String *path = Str_new_from_trusted_utf8(buffer, length);
    Snapshot_Set_Path(snapshot, path);

    buffer += length;

    while(buffer != end)
    {
        size_t length = *((size_t *) buffer);
        buffer += sizeof(size_t);

        String *entry = Str_new_from_trusted_utf8(buffer, length);
        Snapshot_Add_Entry(snapshot, entry);

        buffer += length;
    }

    return snapshot;
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


static void base_commit(CommitRoutineContext *context)
{
    Indexer_Commit(context->lucy->indexer);

    context->binarySnapshot = snapshot_export(Indexer_Get_Snapshot(context->lucy->indexer));

    safeDecref(context->lucy->indexer);
}


bytea *lucy_commit(Lucy *lucy)
{
    CommitRoutineContext context;
    context.lucy = lucy;
    context.binarySnapshot = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_commit, (void *) &context);

    if(error != NULL)
    {
        safeNothrowDecref(context.lucy->indexer);
        throwError(error);
    }

    return context.binarySnapshot;
}


void lucy_rollback(Lucy *lucy)
{
    safeNothrowDecref(lucy->indexer);
}


static void base_set_snapshot(SetSnapshotRoutineContext *context)
{
    if(context->lucy->searcher != NULL)
        safeDecref(context->lucy->searcher);

    context->snapshot = snapshot_import(context->binarySnapshot);
    context->lucy->searcher = IxSearcher_new((Obj *) (context->lucy->folder), context->snapshot);

    safeDecref(context->snapshot);
}


void lucy_set_snapshot(Lucy *lucy, bytea *binarySnapshot)
{
    SetSnapshotRoutineContext context;
    context.lucy = lucy;
    context.binarySnapshot = binarySnapshot;
    context.snapshot = NULL;

    Err *error = Err_trap((Err_Attempt_t) &base_set_snapshot, (void *) &context);

    if(error != NULL)
    {
        safeNothrowDecref(context.lucy->searcher);
        safeNothrowDecref(context.snapshot);

        throwError(error);
    }
}


static void base_search(SearchRoutineContext *context)
{
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
