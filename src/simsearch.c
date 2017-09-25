#include <postgres.h>
#include <access/htup_details.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <fmgr.h>
#include <funcapi.h>
#include <utils/array.h>
#include <utils/builtins.h>
#include <utils/memutils.h>
#include "bitset.h"
#include "heap.h"
#include "java.h"
#include "pgchem.h"
#include "simsearch.h"


typedef struct
{
    int32_t topN;
    float4 cutoff;

    BitSet fp;
    int queryBitCount;

    float4 bound;
    int32_t lowBucketNum;
    int32_t highBucketNum;
    int32_t currBucketNum;
    int32_t foundResults;

    SPITupleTable *table;
    int tableRowCount;
    int tableRowPosition;

    Heap heap;

#if SHOW_STATS
    struct timeval begin;
#endif

} SimilaritySearchData;


static bool initialised = false;
static int moleculeCount;
static SPIPlanPtr mainQueryPlan;
static MemoryContext mcxt = NULL;
static TupleDesc tupdesc = NULL;


void simsearch_module_init(void)
{
    mcxt = AllocSetContextCreate(TopMemoryContext, "simsearch memory context", ALLOCSET_DEFAULT_SIZES);


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "simsearch module: SPI_connect() failed");


    /* load count of molecule records */
    if(unlikely(SPI_execute("select count(*) from " SIMILARITY_FINGERPRINT_TABLE, true, FETCH_ALL) != SPI_OK_SELECT))
        elog(ERROR, "simsearch module: SPI_execute() failed");

    if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
        elog(ERROR, "simsearch module: SPI_execute() failed");

    char isNullFlag;
    moleculeCount = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "simsearch module: SPI_getbinval() failed");

    SPI_freetuptable(SPI_tuptable);


    /* prepare query plan */
    mainQueryPlan = SPI_prepare("select id, fp from " SIMILARITY_FINGERPRINT_TABLE " where bit_count = $1", 1, (Oid[]) { INT4OID });

    if(unlikely(mainQueryPlan == NULL))
        elog(ERROR, "simsearch module: SPI_prepare_cursor() failed");

    if(unlikely(SPI_keepplan(mainQueryPlan) == SPI_ERROR_ARGUMENT))
        elog(ERROR, "simsearch module: SPI_keepplan() failed");


    SPI_finish();
    initialised = true;


    /* create tuple description */
    PG_MEMCONTEXT_BEGIN(mcxt);
    tupdesc = CreateTemplateTupleDesc(2, false);
    TupleDescInitEntry(tupdesc, (AttrNumber) 1, "id", INT4OID, -1, 0);
    TupleDescInitEntry(tupdesc, (AttrNumber) 2, "score", FLOAT4OID, -1, 0);
    tupdesc = BlessTupleDesc(tupdesc);
    PG_MEMCONTEXT_END();
}


void simsearch_module_finish(void)
{
    initialised = false;

    if(mcxt != NULL)
        MemoryContextDelete(mcxt);
}


PG_FUNCTION_INFO_V1(orchem_similarity_search);
Datum orchem_similarity_search(PG_FUNCTION_ARGS)
{
    if(SRF_IS_FIRSTCALL())
    {
#if SHOW_STATS
        struct timeval begin;
        gettimeofday(&begin, NULL);
#endif

        if(unlikely(!initialised))
            elog(ERROR, "simsearch module is not properly initialized");

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        text *type = PG_GETARG_TEXT_P(1);
        float4 cutoff = PG_GETARG_FLOAT4(2);
        int32_t topN = PG_GETARG_INT32(3);
        char *typeStr = text_to_cstring(type);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();
        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);

        SimilaritySearchData *info = palloc(sizeof(SimilaritySearchData));
        funcctx->user_fctx = info;

        info->cutoff = cutoff;
        info->topN = topN;

        uint64_t *words;
        int length =  java_parse_similarity_query(&words, VARDATA(query), VARSIZE(query) - VARHDRSZ, typeStr);

        PG_FREE_IF_COPY(query, 0);
        PG_FREE_IF_COPY(type, 1);
        pfree(typeStr);

        if(likely(length >= 0))
        {
            bitset_init(&info->fp, words, length);
            heap_init(&info->heap, moleculeCount);

            info->queryBitCount = bitset_cardinality(&info->fp);
            info->lowBucketNum = info->queryBitCount - 1;
            info->highBucketNum = info->queryBitCount + 1;
            info->currBucketNum = info->queryBitCount;
            info->bound = 1.0f;
            info->foundResults = 0;

            info->table = NULL;
            info->tableRowCount = -1;
            info->tableRowPosition = -1;
        }
        else
        {
            // fake values to stop searching
            info->topN = INT32_MAX;
            info->foundResults = INT32_MAX;
        }

#if SHOW_STATS
        info->begin = begin;
#endif

        PG_MEMCONTEXT_END();
    }


    bool connected = false;

    FuncCallContext *funcctx = SRF_PERCALL_SETUP();
    SimilaritySearchData *info = funcctx->user_fctx;
    Heap *heap = &info->heap;

    HeapItem result;
    bool isNull = true;


    if(likely(info->topN <= 0 || info->topN != info->foundResults))
    {
        while(true)
        {
            if(heap_size(heap) > 0 && heap_head(heap).score >= info->bound)
            {
                result = heap_head(heap);
                heap_remove(heap);
                isNull = false;
                break;
            }

            if(info->bound < info->cutoff)
                break;

            if(unlikely(info->table == NULL))
            {
                if(unlikely(!connected && SPI_connect() != SPI_OK_CONNECT))
                     elog(ERROR, "simsearch module: SPI_connect() failed");

                connected = true;

                Datum values[] = { Int32GetDatum(info->currBucketNum)};

                if(unlikely(SPI_execute_plan(mainQueryPlan, values, NULL, true, 0) != SPI_OK_SELECT))
                    elog(ERROR, "simsearch module: SPI_execute_plan() failed");

                if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
                    elog(ERROR, "simsearch module: SPI_execute_plan() failed");

                info->table = SPI_tuptable;
                info->tableRowCount = SPI_processed;
                info->tableRowPosition = 0;

                MemoryContextSetParent(SPI_tuptable->tuptabcxt, funcctx->multi_call_memory_ctx);
            }


            if(unlikely(info->tableRowPosition == info->tableRowCount))
            {
                MemoryContextDelete(info->table->tuptabcxt);
                info->table = NULL;

                float up = info->queryBitCount / info->highBucketNum;
                float down = info->lowBucketNum / info->queryBitCount;

                if(up > down)
                {
                    info->currBucketNum = info->highBucketNum;
                    info->highBucketNum++;
                }
                else
                {
                    info->currBucketNum = info->lowBucketNum;
                    info->lowBucketNum--;
                }

                if(info->lowBucketNum < 1 && info->highBucketNum > EXTFP_SIZE)
                    break;

                if(info->currBucketNum < info->queryBitCount)
                    info->bound = info->currBucketNum / (float4) info->queryBitCount;
                else
                    info->bound = info->queryBitCount / (float4) info->currBucketNum;

                continue;
            }


            TupleDesc tupdesc = info->table->tupdesc;
            HeapTuple tuple = info->table->vals[info->tableRowPosition++];
            char isNullFlag;


            Datum id = SPI_getbinval(tuple, tupdesc, 1, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "simsearch module: SPI_getbinval() failed");


            Datum fpDatum = SPI_getbinval(tuple, tupdesc, 2, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "simsearch module: SPI_getbinval() failed");


            BitSet fp;
            ArrayType *fpArray = DatumGetArrayTypeP(fpDatum);
            bitset_init(&fp, (uint64_t *) ARR_DATA_PTR(fpArray), ARR_DIMS(fpArray)[0]);
            int targetBitCount = bitset_cardinality(&fp);

            int bitsInCommon = bitset_and_cardinality(&info->fp, &fp);
            float4 score = bitsInCommon / (float4) (info->queryBitCount + targetBitCount - bitsInCommon);

            if(score == info->bound)
            {
                result.id = id;
                result.score = score;
                isNull = false;
                break;
            }

            if(score >= info->cutoff)
                heap_add(&info->heap, (HeapItem) {.id = DatumGetInt32(id), .score = score});
        }
    }

    if(connected)
        SPI_finish();

    if(unlikely(isNull))
    {
#if SHOW_STATS
        struct timeval begin = ((SimilaritySearchData *) funcctx->user_fctx)->begin;
        struct timeval end;
        gettimeofday(&end, NULL);
        int32_t time_spent = ((int64_t) end.tv_sec - (int64_t) begin.tv_sec) * 1000000 + ((int64_t) end.tv_usec - (int64_t) begin.tv_usec);
        elog(NOTICE, "stat: %i %i.%i ms", info->resultCount, time_spent / 1000, time_spent % 1000);
#endif

        SRF_RETURN_DONE(funcctx);
    }
    else
    {
        info->foundResults++;

        char isnull[2] = {0, 0};
        Datum values[2] = {Int32GetDatum(result.id), Float4GetDatum(result.score)};
        HeapTuple tuple = heap_form_tuple(tupdesc, values, isnull);

        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(tuple));
    }
}
