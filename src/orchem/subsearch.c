#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <funcapi.h>
#include <utils/memutils.h>
#include <sys/mman.h>
#include <unistd.h>
#include "bitset.h"
#include "common.h"
#include "isomorphism.h"
#include "molecule.h"
#include "sachem.h"
#include "subsearch.h"
#include "measurement.h"
#include "java/orchem.h"


#define SHOW_STATS              0
#define FETCH_SIZE              10000


typedef struct
{
    int32_t topN;
    bool extended;
    GraphMode graphMode;
    ChargeMode chargeMode;
    IsotopeMode isotopeMode;
    StereoMode stereoMode;
    int32_t vf2_timeout;

    BitSet resultMask;
    int32_t foundResults;

    OrchemSubstructureQueryData *queryData;
    int queryDataCount;
    int queryDataPosition;

    BitSet candidates;
    int candidatePosition;

    SPITupleTable *table;
    int tableRowCount;
    int tableRowPosition;

    Molecule queryMolecule;
    VF2State vf2state;

    ArrayType *arrayBuffer;
    MemoryContext isomorphismContext;
    MemoryContext targetContext;

#if SHOW_STATS
    int candidateCount;
    struct timeval begin;
    int64_t prepareTime;
    int64_t indexTime;
    int64_t matchTime;
#endif

} SubstructureSearchData;


static bool javaInitialized = false;
static int indexId = -1;
static uint64_t *base = MAP_FAILED;
static size_t length;
static int moleculeCount;
static BitSet bitmap[FP_SIZE];
static SPIPlanPtr indexQueryPlan = NULL;
static SPIPlanPtr mainQueryPlan = NULL;

#if USE_COUNT_FINGERPRINT
static int16 (*counts)[COUNTS_SIZE];
#endif

#if SHOW_STATS
static int dbSize;
#endif


static void orchem_subsearch_init(void)
{
    if(unlikely(javaInitialized == false))
    {
        java_orchem_init();
        javaInitialized = true;
    }


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);


    /* prepare index query plan */
    if(unlikely(indexQueryPlan == NULL))
    {
        SPIPlanPtr plan = SPI_prepare("select id, path from " INDEX_TABLE, 0, NULL);

        if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        indexQueryPlan = plan;
    }


    /* main query plan */
    if(unlikely(mainQueryPlan == NULL))
    {
        SPIPlanPtr plan = SPI_prepare("select id, seqid, molecule from " MOLECULES_TABLE " where seqid = any($1)", 1, (Oid[]) { INT4ARRAYOID });

        if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        mainQueryPlan = plan;
    }


    /* get index information */
    if(unlikely(SPI_execute_plan(indexQueryPlan, NULL, NULL, true, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);


    char isNullFlag;

    Datum id = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);

    if(unlikely(DatumGetInt32(id) != indexId))
    {
        if(likely(base != MAP_FAILED))
        {
            if(unlikely(munmap(base, length) < 0))
                elog(ERROR, "%s: munmap() failed", __func__);

            base = MAP_FAILED;
        }


        Datum path = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        int fd = -1;
        void *address = MAP_FAILED;
        struct stat st;

        PG_TRY();
        {
            if(unlikely((fd = open(TextDatumGetCString(path), O_RDONLY, 0)) < 0))
                elog(ERROR, "%s: open() failed", __func__);

            if(fstat(fd, &st) < 0)
                elog(ERROR, "%s: fstat() failed", __func__);

            if(unlikely((address = mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED))
                elog(ERROR, "%s: mmap() failed", __func__);

            if(unlikely(close(fd) < 0))
                elog(ERROR, "%s: close() failed", __func__);
        }
        PG_CATCH();
        {
            if(address != MAP_FAILED)
                munmap(address, st.st_size);

            if(fd != -1)
                close(fd);

            PG_RE_THROW();
        }
        PG_END_TRY();

        base = address;
        length = st.st_size;
        moleculeCount = *(base + FP_SIZE);

        for(int i = 0; i < FP_SIZE; i++)
        {
            uint64_t *address = base + base[i];
            bitset_init(bitmap + i, address + 1, *address);
        }

#if USE_COUNT_FINGERPRINT
        counts = (int16 (*)[COUNTS_SIZE]) (base + FP_SIZE + 1);
#endif

#if SHOW_STATS
        if(unlikely(SPI_execute("select count(*) from " MOLECULES_TABLE, true, FETCH_ALL) != SPI_OK_SELECT))
            elog(ERROR, "%s: SPI_execute() failed", __func__);

        if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
            elog(ERROR, "%s: SPI_execute() failed", __func__);

        dbSize = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        SPI_freetuptable(SPI_tuptable);
#endif

        indexId = DatumGetInt32(id);
    }

    SPI_finish();
}


PG_FUNCTION_INFO_V1(orchem_substructure_search);
Datum orchem_substructure_search(PG_FUNCTION_ARGS)
{
    if(SRF_IS_FIRSTCALL())
    {
#if SHOW_STATS
        struct timeval begin = time_get();
#endif
        orchem_subsearch_init();

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        int32_t type = PG_GETARG_INT32(1);
        int32_t topN = PG_GETARG_INT32(2);
        GraphMode graphMode = PG_GETARG_INT32(3);
        ChargeMode chargeMode = PG_GETARG_INT32(4);
        IsotopeMode isotopeMode = PG_GETARG_INT32(5);
        StereoMode stereoMode = PG_GETARG_INT32(6);
        TautomerMode tautomerMode = PG_GETARG_INT32(7);
        int32_t vf2_timeout = PG_GETARG_INT32(8);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();
        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);

        SubstructureSearchData *info = (SubstructureSearchData *) palloc(sizeof(SubstructureSearchData));
        funcctx->user_fctx = info;

        info->topN = topN;
        info->graphMode = graphMode;
        info->chargeMode = chargeMode;
        info->isotopeMode = isotopeMode;
        info->stereoMode = stereoMode;
        info->vf2_timeout = vf2_timeout;

#if SHOW_STATS
        struct timeval java_begin = time_get();
#endif
        info->queryDataCount = java_orchem_parse_substructure_query(&info->queryData, VARDATA(query), VARSIZE(query) - VARHDRSZ,
                type, graphMode == GRAPH_EXACT, tautomerMode == TAUTOMER_INCHI);
#if SHOW_STATS
        struct timeval java_end = time_get();
#endif

        PG_FREE_IF_COPY(query, 0);

        info->queryDataPosition = -1;
        info->candidatePosition = -1;
        info->table = NULL;
        info->tableRowCount = -1;
        info->tableRowPosition = -1;
        info->foundResults = 0;

        bitset_init_alloc(&info->candidates, moleculeCount);
        bitset_init_setted(&info->resultMask, moleculeCount);

        info->isomorphismContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx,
                "subsearch isomorphism context", ALLOCSET_DEFAULT_SIZES);
        info->targetContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx,
                "subsearch target context", ALLOCSET_DEFAULT_SIZES);

        info->arrayBuffer = (ArrayType *) palloc(FETCH_SIZE * sizeof(int32) + ARR_OVERHEAD_NONULLS(1));
        info->arrayBuffer->ndim = 1;
        info->arrayBuffer->dataoffset = 0;
        info->arrayBuffer->elemtype = INT4OID;

#if SHOW_STATS
        info->begin = begin;
        info->candidateCount = 0;
        info->prepareTime = time_spent(java_begin, java_end);
        info->indexTime = 0;
        info->matchTime = 0;
#endif

        PG_MEMCONTEXT_END();
    }


    bool connected = false;

    FuncCallContext *funcctx = SRF_PERCALL_SETUP();
    SubstructureSearchData *info = funcctx->user_fctx;

    Datum result;
    bool isNull = true;

    if(likely(info->topN <= 0 || info->topN != info->foundResults))
    {
        while(true)
        {
            if(unlikely(info->tableRowPosition == info->tableRowCount))
            {
                if(info->table != NULL)
                {
                    MemoryContextDelete(info->table->tuptabcxt);
                    info->table = NULL;
                }

                if(info->candidatePosition < 0)
                {
                    info->queryDataPosition++;

                    if(unlikely(info->queryDataPosition == info->queryDataCount))
                        break;

#if SHOW_STATS
                    struct timeval search_begin = time_get();
#endif
                    bitset_copy(&info->candidates, &info->resultMask);

                    for(int i = 0; i < info->queryData[info->queryDataPosition].fpLength; i++)
                    {
                        int idx = Int16GetDatum(info->queryData[info->queryDataPosition].fp[i]);
                        bitset_merge(&info->candidates, bitmap + idx);
                    }
#if SHOW_STATS
                    struct timeval search_end = time_get();
                    info->indexTime += time_spent(search_begin, search_end);
#endif

                    OrchemSubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                    MemoryContextReset(info->isomorphismContext);
                    PG_MEMCONTEXT_BEGIN(info->isomorphismContext);

#if SHOW_STATS
                    struct timeval vf2init_begin = time_get();
#endif
                    info->extended = molecule_is_extended_search_needed(data->molecule, info->chargeMode, info->isotopeMode);
                    molecule_init(&info->queryMolecule, data->molecule, data->restH, info->extended,
                            info->chargeMode != CHARGE_IGNORE, info->isotopeMode != ISOTOPE_IGNORE, info->stereoMode != STEREO_IGNORE);
                    vf2state_init(&info->vf2state, &info->queryMolecule, info->graphMode, info->chargeMode,
                            info->isotopeMode, info->stereoMode);
#if SHOW_STATS
                    struct timeval vf2init_end = time_get();
                    info->prepareTime += time_spent(vf2init_begin, vf2init_end);
#endif

                    PG_MEMCONTEXT_END();

                    info->candidatePosition = bitset_next_set_bit(&info->candidates, 0);
                }


                int32 *arrayData = (int32 *) ARR_DATA_PTR(info->arrayBuffer);
                OrchemSubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                size_t count = 0;

#if SHOW_STATS
                struct timeval get_begin = time_get();
#endif
                while(count < FETCH_SIZE && info->candidatePosition >= 0)
                {
#if USE_COUNT_FINGERPRINT
                    bool isValid = true;

                    for(int j = 0; j < COUNTS_SIZE; j++)
                    {
                        if(counts[info->candidatePosition][j] < data->counts[j])
                        {
                            isValid = false;
                            break;
                        }
                    }

                    if(isValid)
#endif
                        arrayData[count++] = info->candidatePosition;

                    info->candidatePosition = bitset_next_set_bit(&info->candidates, info->candidatePosition + 1);
                }
#if SHOW_STATS
                struct timeval get_end = time_get();
                info->indexTime += time_spent(get_begin, get_end);
#endif

                if(unlikely(count == 0))
                    continue;

#if SHOW_STATS
                struct timeval load_begin = time_get();
#endif
                *(ARR_DIMS(info->arrayBuffer)) = count;
                SET_VARSIZE(info->arrayBuffer, count * sizeof(int32) + ARR_OVERHEAD_NONULLS(1));

                Datum values[] = { PointerGetDatum(info->arrayBuffer)};


                if(unlikely(!connected && SPI_connect() != SPI_OK_CONNECT))
                     elog(ERROR, "%s: SPI_connect() failed", __func__);

                connected = true;


                if(unlikely(SPI_execute_plan(mainQueryPlan, values, NULL, true, 0) != SPI_OK_SELECT))
                    elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

                if(unlikely(SPI_processed != count || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 3))
                    elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

                info->table = SPI_tuptable;
                info->tableRowCount = SPI_processed;
                info->tableRowPosition = 0;

                MemoryContextSetParent(SPI_tuptable->tuptabcxt, funcctx->multi_call_memory_ctx);
#if SHOW_STATS
                struct timeval load_end = time_get();
                info->matchTime += time_spent(load_begin, load_end);
#endif
            }

#if SHOW_STATS
            struct timeval match_begin = time_get();
#endif

            TupleDesc tupdesc = info->table->tupdesc;
            HeapTuple tuple = info->table->vals[info->tableRowPosition++];
            char isNullFlag;


            Datum id = SPI_getbinval(tuple, tupdesc, 1, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);


            Datum seqid = SPI_getbinval(tuple, tupdesc, 2, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);


            Datum molecule = SPI_getbinval(tuple, tupdesc, 3, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);

#if SHOW_STATS
            info->candidateCount++;
#endif

            bytea *moleculeData = DatumGetByteaP(molecule);
            bool match;

            PG_MEMCONTEXT_BEGIN(info->targetContext);
            Molecule target;
            molecule_init(&target, (uint8_t *) VARDATA(moleculeData), NULL, info->extended, info->chargeMode != CHARGE_IGNORE,
                    info->isotopeMode != ISOTOPE_IGNORE, info->stereoMode != STEREO_IGNORE);
            match = vf2state_match(&info->vf2state, &target, DatumGetInt32(id), info->vf2_timeout);
            PG_MEMCONTEXT_END();
            MemoryContextReset(info->targetContext);

#if SHOW_STATS
            struct timeval match_end = time_get();
            info->matchTime += time_spent(match_begin, match_end);
#endif

            if(match)
            {
                bitset_unset(&info->resultMask, DatumGetInt32(seqid));
                info->foundResults++;
                result = id;
                isNull = false;
                break;
            }
        }
    }

    if(connected)
        SPI_finish();

    if(unlikely(isNull))
    {
#if SHOW_STATS
        struct timeval end = time_get();
        int64_t spentTime = time_spent(info->begin, end);
        int64_t sumTime = info->prepareTime + info->indexTime + info->matchTime;
        double scale = spentTime / 100.0;

        double specificity = 100.0 * (dbSize - info->candidateCount) / (dbSize - info->foundResults);
        double precision = 100.0 * info->foundResults / info->candidateCount;

        if(info->candidateCount == 0 && info->foundResults == 0)
            precision = 100.0;

        elog(NOTICE, "stat:   candidates      results  specificity    precision");
        elog(NOTICE, "stat: %12i %12i %11.3f%% %11.3f%%", info->candidateCount, info->foundResults, specificity, precision);
        elog(NOTICE, "time:  fingerprint        index        match          sum        total");
        elog(NOTICE, "time: %10.3fms %10.3fms %10.3fms %10.3fms %10.3fms", time_to_ms(info->prepareTime),
                time_to_ms(info->indexTime), time_to_ms(info->matchTime), time_to_ms(sumTime), time_to_ms(spentTime));
        elog(NOTICE, "time: %11.2f%% %11.2f%% %11.2f%% %11.2f%% ", info->prepareTime / scale,
                info->indexTime / scale, info->matchTime / scale, sumTime / scale);
#endif

        SRF_RETURN_DONE(funcctx);
    }
    else
    {
        SRF_RETURN_NEXT(funcctx, result);
    }
}
