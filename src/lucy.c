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


void lucy_init(Lucy *lucy, const char *index_dir)
{
    lucy_bootstrap_parcel();

    lucy->folder = Str_newf(index_dir);
    lucy->idF = Str_newf("id");
    lucy->fpF = Str_newf("fp");

    lucy->schema = Schema_new();
    StringType *stype = StringType_new();
    StringType_Set_Indexed(stype, false);
    Schema_Spec_Field(lucy->schema, lucy->idF, (FieldType*) stype);
    DECREF(stype);

    String* re = Str_newf ("[a-zA-Z0-9+/]+");
    RegexTokenizer *rt = RegexTokenizer_new(re);
    FullTextType *fttype = FullTextType_new((Analyzer*) rt);
    FullTextType_Set_Indexed(fttype, true);
    FullTextType_Set_Highlightable(fttype, false);
    FullTextType_Set_Stored(fttype, false);
    Schema_Spec_Field(lucy->schema, lucy->fpF, (FieldType*) fttype);
    DECREF(rt);
    DECREF(re);
    DECREF(fttype);

    String *boolop = Str_newf("AND");
    lucy->qparser = QParser_new(lucy->schema, NULL, boolop, NULL);
    DECREF(boolop);
}


void lucy_begin(Lucy *lucy)
{
    lucy->indexer = Indexer_new(lucy->schema, (Obj*) (lucy->folder), NULL, Indexer_CREATE);
}


void lucy_add(Lucy *lucy, int32_t id, const char *fp)
{
    Doc *doc = Doc_new(NULL, 0);

    char buffer[12];
    sprintf(buffer, "%u", id);
    String *idValue = Str_newf(buffer);
    Doc_Store(doc, lucy->idF, (Obj*) idValue);
    DECREF(idValue);

    String *fpValue = Str_newf(fp);
    Doc_Store(doc, lucy->fpF, (Obj*) fpValue);
    DECREF(fpValue);

    Indexer_Add_Doc(lucy->indexer, doc, 1.0);
    DECREF(doc);
}


void lucy_delete(Lucy *lucy, int32_t id)
{
    char buffer[12];
    sprintf(buffer, "%u", id);
    String *idValue = Str_newf(buffer);

    Indexer_Delete_By_Term(lucy->indexer,lucy->idF, (Obj*) idValue);
    DECREF(idValue);
}


void lucy_commit(Lucy *lucy)
{
    Indexer_Commit(lucy->indexer);
    DECREF(lucy->indexer);
}


void lucy_optimize(Lucy *lucy)
{
    Indexer_Optimize(lucy->indexer);
}


Hits *lucy_search(Lucy *lucy, const char *queryFp, int max_results)
{
    lucy->searcher = IxSearcher_new((Obj*) (lucy->folder));

    String *queryStr = Str_newf(queryFp);
    Query *query = QParser_Parse(lucy->qparser, queryStr);
    Hits *hits = IxSearcher_Hits(lucy->searcher, (Obj*) query, 0, max_results, NULL);

    DECREF(query);
    DECREF(queryStr);

    return hits;
}


size_t lucy_get(Lucy *lucy, Hits *hits, int *results, size_t size)
{
    int ret = 0;

    while(size > 0)
    {
        HitDoc *hit = Hits_Next(hits);

        if(hit == NULL)
        {
            DECREF(lucy->searcher);
            break;
        }

        String *id = (String*) HitDoc_Extract(hit, lucy->idF);
        *(results++) = Str_To_I64(id);

        ++ret;
        --size;

        DECREF(id);
        DECREF(hit);
    }

    return ret;
}
