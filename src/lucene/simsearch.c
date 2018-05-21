#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <funcapi.h>
#include <utils/memutils.h>
#include <access/htup_details.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include "common.h"
#include "search.h"
#include "molecule.h"
#include "sachem.h"
#include "lucene.h"
#include "measurement.h"
#include "java/parse.h"
#include "fingerprints/fingerprint.h"


#define SHOW_STATS              0


typedef struct
{
    LuceneSimsearchResult result;

#if SHOW_STATS
    int32_t foundResults;
    struct timeval begin;
#endif

} SimilaritySearchData;


static bool initialized = false;
static TupleDesc tupdesc = NULL;


void lucene_simsearch_init(void)
{
    if(unlikely(initialized == false))
    {
        lucene_search_init();


        /* create tuple description */
        if(unlikely(tupdesc == NULL))
        {
            if(unlikely(SPI_connect() != SPI_OK_CONNECT))
                elog(ERROR, "%s: SPI_connect() failed", __func__);

            TupleDesc desc = NULL;

            PG_MEMCONTEXT_BEGIN(TopMemoryContext);
            PG_TRY();
            {
                desc = CreateTemplateTupleDesc(2, false);
                TupleDescInitEntry(desc, (AttrNumber) 1, "compound", INT4OID, -1, 0);
                TupleDescInitEntry(desc, (AttrNumber) 2, "score", FLOAT4OID, -1, 0);
                desc = BlessTupleDesc(desc);
                tupdesc = desc;
            }
            PG_CATCH();
            {
                if(desc != NULL)
                    FreeTupleDesc(desc);

                PG_RE_THROW();
            }
            PG_END_TRY();
            PG_MEMCONTEXT_END();

            SPI_finish();
        }

        initialized = true;
    }


    /* get snapshot information */
    lucene_search_update_snapshot();
}


PG_FUNCTION_INFO_V1(lucene_similarity_search);
Datum lucene_similarity_search(PG_FUNCTION_ARGS)
{
    if(SRF_IS_FIRSTCALL())
    {
#if SHOW_STATS
        struct timeval begin = time_get();
#endif
        lucene_simsearch_init();

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        int32_t type = PG_GETARG_INT32(1);
        int32_t topN = PG_GETARG_INT32(3);
        float4 cutoff = PG_GETARG_FLOAT4(2);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();

        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);

        SimilaritySearchData *info = (SimilaritySearchData *) palloc(sizeof(SimilaritySearchData));
        funcctx->user_fctx = info;

        SimilarityQueryData queryData;
        java_parse_similarity_query(&queryData, VARDATA(query), VARSIZE(query) - VARHDRSZ, type);

        Molecule molecule;
        molecule_simple_init(&molecule, queryData.molecule);

        IntegerFingerprint fp = integer_similarity_fingerprint_get_query(&molecule);

        info->result = lucene_simsearch_submit(&lucene, fp, topN, cutoff);

        PG_FREE_IF_COPY(query, 0);

#if SHOW_STATS
        info->foundResults = 0;
        info->begin = begin;
#endif

        PG_MEMCONTEXT_END();
    }


    FuncCallContext *funcctx = SRF_PERCALL_SETUP();
    SimilaritySearchData *info = funcctx->user_fctx;

    int32_t id;
    float score;
    bool isNull = true;

    PG_TRY();
    {
        isNull = !lucene_simsearch_get(&lucene, &info->result, &id, &score);
    }
    PG_CATCH();
    {
        lucene_simsearch_fail(&lucene, &info->result);

        PG_RE_THROW();
    }
    PG_END_TRY();

    if(unlikely(isNull))
    {
#if SHOW_STATS
        struct timeval end = time_get();
        int64_t spentTime = time_spent(info->begin, end);
        elog(NOTICE, "stat: %i", info->foundResults);
        elog(NOTICE, "time: %.3fms", time_to_ms(spentTime));
#endif

        SRF_RETURN_DONE(funcctx);
    }
    else
    {
#if SHOW_STATS
        info->foundResults++;
#endif

        char isnull[2] = {0, 0};
        Datum values[2] = {Int32GetDatum(id), Float4GetDatum(score)};
        HeapTuple tuple = heap_form_tuple(tupdesc, values, isnull);

        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(tuple));
    }
}
