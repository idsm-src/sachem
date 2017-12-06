#include <postgres.h>
#include <access/htup.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <fmgr.h>
#include <funcapi.h>
#include <utils/array.h>
#include <utils/builtins.h>
#include <utils/memutils.h>
#include <access/htup_details.h>
#include "bitset.h"
#include "heap.h"
#include "java.h"
#include "sachem.h"
#include "simsearch.h"


#define SHOW_STATS          0
#define EXTFP_SIZE          907
#define FINGERPRINT_TABLE   "orchem_fingerprint"


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
static SPIPlanPtr mainQueryPlan;
static MemoryContext mcxt = NULL;
static TupleDesc tupdesc = NULL;


void simsearch_module_init(void)
{
    PG_TRY();
    {
        /* prepare query plan */
        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
            elog(ERROR, "%s: SPI_connect() failed", __func__);

        mainQueryPlan = SPI_prepare("select id, fp from " FINGERPRINT_TABLE " where bit_count = $1", 1, (Oid[]) { INT4OID });

        if(unlikely(mainQueryPlan == NULL))
            elog(ERROR, "%s: SPI_prepare_cursor() failed", __func__);

        if(unlikely(SPI_keepplan(mainQueryPlan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        SPI_finish();


        /* create tuple description */
        mcxt = AllocSetContextCreate(TopMemoryContext, "simsearch memory context", ALLOCSET_DEFAULT_SIZES);

        PG_MEMCONTEXT_BEGIN(mcxt);
        tupdesc = CreateTemplateTupleDesc(2, false);
        TupleDescInitEntry(tupdesc, (AttrNumber) 1, "compound", INT4OID, -1, 0);
        TupleDescInitEntry(tupdesc, (AttrNumber) 2, "score", FLOAT4OID, -1, 0);
        tupdesc = BlessTupleDesc(tupdesc);
        PG_MEMCONTEXT_END();


        initialised = true;
    }
    PG_CATCH();
    {
        elog(NOTICE, "%s: initialization failed", __func__);
    }
    PG_END_TRY();
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
            elog(ERROR, "%s: simsearch module is not properly initialized", __func__);

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        text *type = PG_GETARG_TEXT_P(1);
        float4 cutoff = PG_GETARG_FLOAT4(2);
        int32_t topN = PG_GETARG_INT32(3);
        char *typeStr = text_to_cstring(type);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();
        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);

        SimilaritySearchData *info = (SimilaritySearchData *) palloc(sizeof(SimilaritySearchData));
        funcctx->user_fctx = info;

        info->cutoff = cutoff;
        info->topN = topN;

        uint64_t *words;
        int length =  java_parse_orchem_similarity_query(&words, VARDATA(query), VARSIZE(query) - VARHDRSZ, typeStr);

        PG_FREE_IF_COPY(query, 0);
        PG_FREE_IF_COPY(type, 1);
        pfree(typeStr);

        if(likely(length >= 0))
        {
            bitset_init(&info->fp, words, length);
            heap_init(&info->heap);

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
                     elog(ERROR, "%s: SPI_connect() failed", __func__);

                connected = true;

                Datum values[] = { Int32GetDatum(info->currBucketNum)};

                if(unlikely(SPI_execute_plan(mainQueryPlan, values, NULL, true, 0) != SPI_OK_SELECT))
                    elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

                if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
                    elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

                info->table = SPI_tuptable;
                info->tableRowCount = SPI_processed;
                info->tableRowPosition = 0;

                MemoryContextSetParent(SPI_tuptable->tuptabcxt, funcctx->multi_call_memory_ctx);
            }


            if(unlikely(info->tableRowPosition == info->tableRowCount))
            {
                MemoryContextDelete(info->table->tuptabcxt);
                info->table = NULL;

                float up = info->queryBitCount / (float) info->highBucketNum;
                float down = info->lowBucketNum / (float) info->queryBitCount;

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
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);


            Datum fpDatum = SPI_getbinval(tuple, tupdesc, 2, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);


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
