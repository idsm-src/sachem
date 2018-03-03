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
#include "java.h"
#include "molecule.h"
#include "sachem.h"
#include "subsearch.h"


#define SHOW_STATS              0
#define FETCH_SIZE              1000


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
#endif

} SubstructureSearchData;


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


static void orchem_subsearch_init(void)
{
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

        indexId = DatumGetInt32(id);
    }

    SPI_finish();
}


PG_FUNCTION_INFO_V1(orchem_substructure_search);
Datum orchem_substructure_search(PG_FUNCTION_ARGS)
{
    if(SRF_IS_FIRSTCALL())
    {
        orchem_subsearch_init();

#if SHOW_STATS
        struct timeval begin;
        gettimeofday(&begin, NULL);
#endif

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

        info->queryDataCount = java_parse_orchem_substructure_query(&info->queryData, VARDATA(query), VARSIZE(query) - VARHDRSZ,
                type, graphMode == GRAPH_EXACT, tautomerMode == TAUTOMER_INCHI);

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


                    bitset_copy(&info->candidates, &info->resultMask);

                    for(int i = 0; i < info->queryData[info->queryDataPosition].fpLength; i++)
                    {
                        int idx = Int16GetDatum(info->queryData[info->queryDataPosition].fp[i]);
                        bitset_merge(&info->candidates, bitmap + idx);
                    }

                    OrchemSubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                    MemoryContextReset(info->isomorphismContext);

                    PG_MEMCONTEXT_BEGIN(info->isomorphismContext);
                    info->extended = molecule_is_extended_search_needed(data->molecule, info->chargeMode, info->isotopeMode);
                    molecule_init(&info->queryMolecule, data->molecule, data->restH, info->extended,
                            info->chargeMode != CHARGE_IGNORE, info->isotopeMode != ISOTOPE_IGNORE, info->stereoMode != STEREO_IGNORE);
                    vf2state_init(&info->vf2state, &info->queryMolecule, info->graphMode, info->chargeMode,
                            info->isotopeMode, info->stereoMode);
                    PG_MEMCONTEXT_END();

                    info->candidatePosition = bitset_next_set_bit(&info->candidates, 0);
                }


                int32 *arrayData = (int32 *) ARR_DATA_PTR(info->arrayBuffer);
                OrchemSubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                int count = 0;

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

                if(unlikely(count == 0))
                    continue;


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
            }

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
            molecule_init(&target, VARDATA(moleculeData), NULL, info->extended, info->chargeMode != CHARGE_IGNORE,
                    info->isotopeMode != ISOTOPE_IGNORE, info->stereoMode != STEREO_IGNORE);
            match = vf2state_match(&info->vf2state, &target, DatumGetInt32(id), info->vf2_timeout);
            PG_MEMCONTEXT_END();
            MemoryContextReset(info->targetContext);

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
        struct timeval begin = ((SubstructureSearchData *) funcctx->user_fctx)->begin;
        struct timeval end;
        gettimeofday(&end, NULL);
        int32_t time_spent = ((int64_t) end.tv_sec - (int64_t) begin.tv_sec) * 1000000 + ((int64_t) end.tv_usec - (int64_t) begin.tv_usec);
        elog(NOTICE, "stat: %i %i %i.%i ms", info->candidateCount, info->foundResults, time_spent / 1000, time_spent % 1000);
#endif

        SRF_RETURN_DONE(funcctx);
    }
    else
    {
        SRF_RETURN_NEXT(funcctx, result);
    }
}
