#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <fmgr.h>
#include <funcapi.h>
#include <utils/array.h>
#include <utils/builtins.h>
#include <utils/memutils.h>
#include "bitset.h"
#include "isomorphism.h"
#include "java.h"
#include "molecule.h"
#include "pgchem.h"
#include "subsearch.h"


typedef struct
{
    int32_t topN;
    bool strictStereo;
    bool exact;

    BitSet resultMask;
    int32_t foundResults;

    QueryData *queryData;
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


static bool initialised = false;
static int moleculeCount;
static BitSet bitmap[FP_SIZE];
static SPIPlanPtr mainQueryPlan;
static MemoryContext mcxt = NULL;

#if USE_COUNT_FINGERPRINT
static int16 *counts[COUNTS_SIZE];
#endif


void subsearch_module_init(void)
{
    mcxt = AllocSetContextCreate(TopMemoryContext, "subsearch memory context", ALLOCSET_DEFAULT_SIZES);

    char isNullFlag;


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "subsearch module: SPI_connect() failed");


    /* load count of molecule records */
    if(unlikely(SPI_execute("select count(*) from " MOLECULES_TABLE, true, FETCH_ALL) != SPI_OK_SELECT))
        elog(ERROR, "subsearch module: SPI_execute() failed");

    if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
        elog(ERROR, "subsearch module: SPI_execute() failed");

    moleculeCount = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "subsearch module: SPI_getbinval() failed");

    SPI_freetuptable(SPI_tuptable);


#if USE_COUNT_FINGERPRINT
    /* load count fingerprint */
    if(unlikely(SPI_execute("select seqid, molTripleBondCount, molSCount, molOCount, molNCount, molFCount, molClCount, molBrCount, "
            "molICount, molCCount, molPCount from " MOLECULE_COUNTS_TABLE " order by seqid", true, FETCH_ALL) != SPI_OK_SELECT))
        elog(ERROR, "subsearch module: SPI_execute() failed");

    if(unlikely(SPI_processed != moleculeCount || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != COUNTS_SIZE + 1))
        elog(ERROR, "subsearch module: SPI_execute() failed");

    PG_MEMCONTEXT_BEGIN(mcxt);

    for(int j = 0; j < COUNTS_SIZE; j++)
        counts[j] = palloc0(moleculeCount * sizeof(int16));

    for(int i = 0; i < SPI_processed; i++)
    {
        HeapTuple tuple = SPI_tuptable->vals[i];

        int32 seqid = DatumGetInt32(SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag));

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "subsearch module: SPI_getbinval() failed");

        for(int j = 0; j < COUNTS_SIZE; j++)
        {
            int16 value = DatumGetInt16(SPI_getbinval(tuple, SPI_tuptable->tupdesc, j + 2, &isNullFlag));

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch module: SPI_getbinval() failed");

            counts[j][seqid] = value;
        }
    }

    PG_MEMCONTEXT_END();

    SPI_freetuptable(SPI_tuptable);
#endif


    /* load bit fingerprint */
    if(unlikely(SPI_execute("select bitmap from " FINGERPRINT_INDEX_TABLE " order by idx", true, FETCH_ALL) != SPI_OK_SELECT))
        elog(ERROR, "subsearch module: SPI_execute() failed");

    if(unlikely(SPI_processed != FP_SIZE || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
        elog(ERROR, "subsearch module: SPI_execute() failed");

    PG_MEMCONTEXT_BEGIN(mcxt);

    for(int i = 0; i < SPI_processed; i++)
    {
        HeapTuple tuple = SPI_tuptable->vals[i];

        bytea *fp = DatumGetByteaP(SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag));

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "subsearch module: SPI_getbinval() failed");

        size_t length = (VARSIZE(fp) - VARHDRSZ);
        BitSet line;
        bitset_init_from_array(bitmap + i, VARDATA(fp), length);
    }

    PG_MEMCONTEXT_END();

    SPI_freetuptable(SPI_tuptable);


    /* prepare query plan */
    mainQueryPlan = SPI_prepare("select id, seqid, atoms, bonds from " MOLECULES_TABLE " where seqid = any($1)", 1, (Oid[]) { INT4ARRAYOID });

    if(unlikely(mainQueryPlan == NULL))
        elog(ERROR, "subsearch module: SPI_prepare_cursor() failed");

    if(unlikely(SPI_keepplan(mainQueryPlan) == SPI_ERROR_ARGUMENT))
        elog(ERROR, "subsearch module: SPI_keepplan() failed");


    SPI_finish();
    initialised = true;
}


void subsearch_module_finish(void)
{
    initialised = false;

    if(mcxt != NULL)
        MemoryContextDelete(mcxt);
}


PG_FUNCTION_INFO_V1(orchem_substructure_search);
Datum orchem_substructure_search(PG_FUNCTION_ARGS)
{
    if(SRF_IS_FIRSTCALL())
    {
#if SHOW_STATS
        struct timeval begin;
        gettimeofday(&begin, NULL);
#endif

        if(unlikely(!initialised))
            elog(ERROR, "subsearch module is not properly initialized");

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        text *type = PG_GETARG_TEXT_P(1);
        int32_t topN = PG_GETARG_INT32(2);
        bool strictStereo = PG_GETARG_BOOL(3);
        bool exact = PG_GETARG_BOOL(4);
        bool tautomers = PG_GETARG_BOOL(5);
        char *typeStr = text_to_cstring(type);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();
        PG_MEMCONTEXT_BEGIN(MemoryContextSwitchTo(funcctx->multi_call_memory_ctx));

        SubstructureSearchData *info = palloc(sizeof(SubstructureSearchData));
        funcctx->user_fctx = info;

        info->topN = topN;
        info->strictStereo = strictStereo;
        info->exact = exact;

        info->queryDataCount = java_parse_query(&info->queryData, VARDATA(query), VARSIZE(query) - VARHDRSZ, typeStr, tautomers);

        PG_FREE_IF_COPY(query, 0);
        PG_FREE_IF_COPY(type, 1);
        pfree(typeStr);

        info->queryDataPosition = -1;
        info->candidatePosition = -1;
        info->table = NULL;
        info->tableRowCount = -1;
        info->tableRowPosition = -1;
        info->foundResults = 0;

        bitset_init_alloc(&info->candidates, moleculeCount);
        bitset_init_setted(&info->resultMask, moleculeCount);

        info->isomorphismContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx, "subsearch isomorphism context", ALLOCSET_DEFAULT_SIZES);
        info->targetContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx, "subsearch target context", ALLOCSET_DEFAULT_SIZES);

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
                    MemoryContextDelete(info->table->tuptabcxt);


                while(info->candidatePosition < 0)
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

                    QueryData *data = &(info->queryData[info->queryDataPosition]);

                    MemoryContextReset(info->isomorphismContext);

                    PG_MEMCONTEXT_BEGIN(info->isomorphismContext);
                    molecule_init(&info->queryMolecule, data->atomLength, data->atoms, data->bondLength, data->bonds, data->restH, !info->exact);
                    vf2state_init(&info->vf2state, &info->queryMolecule, info->strictStereo, info->exact);
                    PG_MEMCONTEXT_END();

                    info->candidatePosition = bitset_next_set_bit(&info->candidates, 0);
                }


                int32 *arrayAata = (int32 *) ARR_DATA_PTR(info->arrayBuffer);
                QueryData *data = &(info->queryData[info->queryDataPosition]);

                int count = 0;

                while(count < FETCH_SIZE && info->candidatePosition >= 0)
                {
#if USE_COUNT_FINGERPRINT
                    bool isValid = true;

                    for(int j = 0; j < COUNTS_SIZE; j++)
                    {
                        if(counts[j][info->candidatePosition] < data->counts[j])
                        {
                            isValid = false;
                            break;
                        }
                    }

                    if(isValid)
#endif
                        arrayAata[count++] = Int32GetDatum(info->candidatePosition);

                    info->candidatePosition = bitset_next_set_bit(&info->candidates, info->candidatePosition + 1);
                }

                if(unlikely(count == 0))
                    break;


                *(ARR_DIMS(info->arrayBuffer)) = count;
                SET_VARSIZE(info->arrayBuffer, count * sizeof(int32) + ARR_OVERHEAD_NONULLS(1));

                Datum values[] = { PointerGetDatum(info->arrayBuffer)};


                if(unlikely(!connected && SPI_connect() != SPI_OK_CONNECT))
                     elog(ERROR, "subsearch module: SPI_connect() failed");

                connected = true;


                if(unlikely(SPI_execute_plan(mainQueryPlan, values, NULL, true, 0) != SPI_OK_SELECT))
                    elog(ERROR, "subsearch module: SPI_execute_plan() failed");

                if(unlikely(SPI_processed != count || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 4))
                    elog(ERROR, "subsearch module: SPI_execute_plan() failed");

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
                elog(ERROR, "subsearch module: SPI_getbinval() failed");


            Datum seqid = SPI_getbinval(tuple, tupdesc, 2, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch module: SPI_getbinval() failed");


            Datum atoms = SPI_getbinval(tuple, tupdesc, 3, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch module: SPI_getbinval() failed");


            Datum bonds = SPI_getbinval(tuple, tupdesc, 4, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch module: SPI_getbinval() failed");

    #if SHOW_STATS
            info->candidateCount++;
    #endif

            bytea *atomsData = DatumGetByteaP(atoms);
            bytea *bondsData = DatumGetByteaP(bonds);

            int atomsize = (VARSIZE(atomsData) - VARHDRSZ) / ATOM_BLOCK_SIZE;
            int bondsize = (VARSIZE(bondsData) - VARHDRSZ) / BOND_BLOCK_SIZE;

            if(atomsize < info->queryMolecule.atomCount)
                continue;

            if(bondsize < info->queryMolecule.bondCount)
                continue;

            MemoryContext old = MemoryContextSwitchTo(info->targetContext);

            Molecule target;
            molecule_init(&target, VARSIZE(atomsData) - VARHDRSZ, VARDATA(atomsData), VARSIZE(bondsData) - VARHDRSZ, VARDATA(bondsData), NULL, false);

            bool match = vf2state_match(&info->vf2state, &target);

            MemoryContextSwitchTo(old);
            MemoryContextReset(info->targetContext);

            if(match)
            {
                bitset_unset(&info->resultMask, seqid);
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
