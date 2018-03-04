#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <funcapi.h>
#include <utils/memutils.h>
#include "bitset.h"
#include "common.h"
#include "isomorphism.h"
#include "molecule.h"
#include "sachem.h"
#include "subsearch.h"
#include "lucene.h"
#include "java/parse.h"
#include "fingerprints/fingerprint.h"


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

    SubstructureQueryData *queryData;
    int queryDataCount;
    int queryDataPosition;

    LuceneResultSet resultSet;

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


static bool javaInitialized = false;
static bool luceneInitialised = false;
static int indexId = -1;
static SPIPlanPtr mainQueryPlan;
static SPIPlanPtr snapshotQueryPlan;
static Lucene lucene;


void lucene_subsearch_init(void)
{
    if(unlikely(javaInitialized == false))
    {
        java_parse_init();
        javaInitialized = true;
    }


    /* prepare lucene */
    if(unlikely(luceneInitialised == false))
    {
        lucene_init(&lucene);
        luceneInitialised = true;
    }


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    /* prepare snapshot query plan */
    if(unlikely(snapshotQueryPlan == NULL))
    {
        SPIPlanPtr plan = SPI_prepare("select id, path from " INDEX_TABLE, 0, NULL);

        if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        snapshotQueryPlan = plan;
    }


    /* prepare query plan */
    if(unlikely(mainQueryPlan == NULL))
    {
        SPIPlanPtr plan = SPI_prepare("select id, molecule from " MOLECULES_TABLE " where id = any($1)", 1, (Oid[]) { INT4ARRAYOID });

        if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        mainQueryPlan = plan;
    }


    /* get snapshot information */
    if(unlikely(SPI_execute_plan(snapshotQueryPlan, NULL, NULL, true, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);


    char isNullFlag;

    Datum id = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);

    if(unlikely(DatumGetInt32(id) != indexId))
    {
        Datum path = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        lucene_set_folder(&lucene, TextDatumGetCString(path));
        indexId = DatumGetInt32(id);
    }

    SPI_finish();
}


PG_FUNCTION_INFO_V1(lucene_substructure_search);
Datum lucene_substructure_search(PG_FUNCTION_ARGS)
{
    bool connected = false;


    if(SRF_IS_FIRSTCALL())
    {
        lucene_subsearch_init();

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


        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
             elog(ERROR, "%s: SPI_connect() failed", __func__);

        connected = true;

        if(unlikely(SPI_execute("select max(id) + 1 from " MOLECULES_TABLE, true, FETCH_ALL) != SPI_OK_SELECT))
            elog(ERROR, "%s: SPI_execute() failed", __func__);

        if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
            elog(ERROR, "%s: SPI_execute() failed", __func__);

        char isNullFlag;
        int64_t moleculeCount = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        SPI_freetuptable(SPI_tuptable);


        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);

        SubstructureSearchData *info = (SubstructureSearchData *) palloc(sizeof(SubstructureSearchData));
        funcctx->user_fctx = info;

        info->topN = topN;
        info->graphMode = graphMode;
        info->chargeMode = chargeMode;
        info->isotopeMode = isotopeMode;
        info->stereoMode = stereoMode;
        info->vf2_timeout = vf2_timeout;

        info->queryDataCount = java_parse_substructure_query(&info->queryData, VARDATA(query), VARSIZE(query) - VARHDRSZ,
                type, graphMode == GRAPH_EXACT, tautomerMode == TAUTOMER_INCHI);

        PG_FREE_IF_COPY(query, 0);

        info->queryDataPosition = -1;
        info->resultSet = NULL_RESULT_SET;
        info->table = NULL;
        info->tableRowCount = -1;
        info->tableRowPosition = -1;
        info->foundResults = 0;

        bitset_init_empty(&info->resultMask, moleculeCount);

        info->isomorphismContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx,
                "subsearch-lucene isomorphism context", ALLOCSET_DEFAULT_SIZES);
        info->targetContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx,
                "subsearch-lucene target context", ALLOCSET_DEFAULT_SIZES);

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


    FuncCallContext *funcctx = SRF_PERCALL_SETUP();
    SubstructureSearchData *info = funcctx->user_fctx;

    Datum result;
    bool isNull = true;

    PG_TRY();
    {
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


                    if(!lucene_is_open(&lucene, &info->resultSet))
                    {
                        info->queryDataPosition++;

                        if(unlikely(info->queryDataPosition == info->queryDataCount))
                            break;

                        SubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                        MemoryContextReset(info->isomorphismContext);

                        PG_MEMCONTEXT_BEGIN(info->isomorphismContext);
                        info->extended = molecule_is_extended_search_needed(data->molecule, info->chargeMode, info->isotopeMode);
                        molecule_init(&info->queryMolecule, data->molecule, data->restH, info->extended,
                                info->chargeMode != CHARGE_IGNORE, info->isotopeMode != ISOTOPE_IGNORE, info->stereoMode != STEREO_IGNORE);
                        vf2state_init(&info->vf2state, &info->queryMolecule, info->graphMode, info->chargeMode, info->isotopeMode,
                                info->stereoMode);

                        IntegerFingerprint fp = integer_fingerprint_get_query(&info->queryMolecule);
                        info->resultSet = lucene_search(&lucene, fp, INT32_MAX);


                        if(fp.data != NULL)
                            pfree(fp.data);

                        PG_MEMCONTEXT_END();
                    }


                    int32 *arrayData = (int32 *) ARR_DATA_PTR(info->arrayBuffer);
                    SubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                    int count = lucene_get(&lucene, &info->resultSet, arrayData, FETCH_SIZE);


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

                    if(unlikely(SPI_processed != count || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
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


                Datum molecule = SPI_getbinval(tuple, tupdesc, 2, &isNullFlag);

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
                    bitset_set(&info->resultMask, DatumGetInt32(id));
                    info->foundResults++;
                    result = id;
                    isNull = false;
                    break;
                }
            }
        }
    }
    PG_CATCH();
    {
        lucene_fail(&lucene, &info->resultSet);

        PG_RE_THROW();
    }
    PG_END_TRY();

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
