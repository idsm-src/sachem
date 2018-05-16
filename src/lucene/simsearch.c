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
static bool javaInitialized = false;
static bool luceneInitialised = false;
static int indexId = -1;
static SPIPlanPtr snapshotQueryPlan;
static Lucene lucene;
static TupleDesc tupdesc = NULL;


void lucene_simsearch_init(void)
{
    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);


    if(unlikely(initialized == false))
    {
        if(unlikely(javaInitialized == false))
        {
            java_parse_init();
            javaInitialized = true;
        }


        /* create tuple description */
        if(unlikely(tupdesc == NULL))
        {
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
        }


        /* prepare lucene */
        if(unlikely(luceneInitialised == false))
        {
            lucene_init(&lucene);
            luceneInitialised = true;
        }


        /* prepare snapshot query plan */
        if(unlikely(snapshotQueryPlan == NULL))
        {
            SPIPlanPtr plan = SPI_prepare("select id from " INDEX_TABLE, 0, NULL);

            if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
                elog(ERROR, "%s: SPI_keepplan() failed", __func__);

            snapshotQueryPlan = plan;
        }


        initialized = true;
    }


    /* get snapshot information */
    if(unlikely(SPI_execute_plan(snapshotQueryPlan, NULL, NULL, true, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    char isNullFlag;
    int32_t dbIndexNumber = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);


    if(unlikely(dbIndexNumber != indexId))
    {
        char *path = get_index_path(LUCENE_INDEX_PREFIX, LUCENE_INDEX_SUFFIX, dbIndexNumber);
        lucene_set_folder(&lucene, path);
        indexId = dbIndexNumber;
    }

    SPI_finish();
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

        IntegerFingerprint fp = integer_fingerprint_get(&molecule);

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
