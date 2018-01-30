#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <fmgr.h>
#include <funcapi.h>
#include <utils/array.h>
#include <utils/builtins.h>
#include <utils/memutils.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include "bitset.h"
#include "isomorphism.h"
#include "java.h"
#include "molecule.h"
#include "sachem.h"
#include "subsearch.h"
#include "lucy.h"
#include "fingerprints/fingerprint.h"


#define SHOW_STATS                0
#define FETCH_SIZE                1000
#define SYNC_FETCH_SIZE           100000
#define COMPOUNDS_TABLE           "compounds"
#define MOLECULES_TABLE           "sachem_molecules"
#define MOLECULE_ERRORS_TABLE     "sachem_molecule_errors"
#define AUDIT_TABLE               "sachem_compound_audit"


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

    Hits *lucySearch;

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


typedef struct
{
    pg_atomic_uint32 possition;
    uint32_t count;
    Fingerprint *result;
    LucyLoaderData *data;
} ThreadData;


static bool initialised = false;
static SPIPlanPtr mainQueryPlan;
static Lucy lucy;


void subsearch_lucy_module_init(void)
{
    PG_TRY();
    {
        /* prepare lucy */
        lucy_init(&lucy, getFilePath("lucy"));


        /* prepare query plan */
        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
            elog(ERROR, "%s: SPI_connect() failed", __func__);

        mainQueryPlan = SPI_prepare("select id, molecule from " MOLECULES_TABLE " where id = any($1)", 1, (Oid[]) { INT4ARRAYOID });

        if(unlikely(SPI_keepplan(mainQueryPlan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);


        SPI_finish();
        initialised = true;
    }
    PG_CATCH();
    {
        elog(NOTICE, "%s: initialization failed", __func__);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


void subsearch_lucy_module_finish(void)
{
    initialised = false;
}


PG_FUNCTION_INFO_V1(lucy_substructure_search);
Datum lucy_substructure_search(PG_FUNCTION_ARGS)
{
    bool connected = false;


    if(SRF_IS_FIRSTCALL())
    {
#if SHOW_STATS
        struct timeval begin;
        gettimeofday(&begin, NULL);
#endif

        if(unlikely(!initialised))
            elog(ERROR, "%s: subsearch-lucy module is not properly initialized", __func__);

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
        info->lucySearch = NULL;
        info->table = NULL;
        info->tableRowCount = -1;
        info->tableRowPosition = -1;
        info->foundResults = 0;

        bitset_init_empty(&info->resultMask, moleculeCount);

        info->isomorphismContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx,
                "subsearch-lucy isomorphism context", ALLOCSET_DEFAULT_SIZES);
        info->targetContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx,
                "subsearch-lucy target context", ALLOCSET_DEFAULT_SIZES);

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


                if(info->lucySearch == NULL)
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

                    Fingerprint fp = fingerprint_get_query(&info->queryMolecule, &palloc);

                    if(!fingerprint_is_valid(fp))
                        elog(ERROR, "fingerprint cannot be generated");

                    info->lucySearch = lucy_search(&lucy, fp, INT32_MAX);

                    if(fp.data != NULL)
                        pfree(fp.data);

                    PG_MEMCONTEXT_END();
                }


                int32 *arrayData = (int32 *) ARR_DATA_PTR(info->arrayBuffer);
                SubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                int count = lucy_get(&lucy, info->lucySearch, arrayData, FETCH_SIZE);


                if(unlikely(count == 0))
                {
                    info->lucySearch = NULL;
                    continue;
                }


                *(ARR_DIMS(info->arrayBuffer)) = count;
                SET_VARSIZE(info->arrayBuffer, count * sizeof(int32) + ARR_OVERHEAD_NONULLS(1));

                Datum values[] = { PointerGetDatum(info->arrayBuffer)};


                if(unlikely(!connected && SPI_connect() != SPI_OK_CONNECT))
                     elog(ERROR, "%s: SPI_connect() failed", __func__);

                connected = true;


                if(unlikely(SPI_execute_plan(mainQueryPlan, values, NULL, true, 0) != SPI_OK_SELECT))
                    elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

                if(unlikely(/*SPI_processed != count ||*/ SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
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
            match = vf2state_match(&info->vf2state, &target, info->vf2_timeout);
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


void *lucy_substructure_process_data(void *dataPtr)
{
    ThreadData* data = (ThreadData*) dataPtr;
    Molecule molecule;

    while(true)
    {
        uint32_t i = pg_atomic_fetch_add_u32(&data->possition, 1);

        if(i >= data->count)
            break;

        if(data->data[i].molecule)
        {
            if(molecule_simple_init(&molecule, VARDATA(data->data[i].molecule), malloc))
                data->result[i] = fingerprint_get(&molecule, &malloc);

            molecule_simple_free(&molecule, free);
        }
        else
        {
            data->result[i] = (Fingerprint) {.size = -1, .data = NULL};
        }
    }

    pthread_exit((void*) 0);
}


PG_FUNCTION_INFO_V1(lucy_sync_data);
Datum lucy_sync_data(PG_FUNCTION_ARGS)
{
    bool verbose = PG_GETARG_BOOL(0);

    MemoryContext memoryContext = CurrentMemoryContext;
    createBasePath();


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    char isNullFlag;

    /*
     * perepare reqired thread data
     */
    int countOfThread = sysconf(_SC_NPROCESSORS_ONLN);

    pthread_attr_t attr;
    pthread_attr_init(&attr);


    lucy_begin(&lucy);

    /*
     * delete unnecessary data
     */

    Portal auditCursor = SPI_cursor_open_with_args(NULL, "select id from " AUDIT_TABLE " where not stored",
            0, NULL, NULL, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

    while(true)
    {
        SPI_cursor_fetch(auditCursor, true, SYNC_FETCH_SIZE);

        if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
            elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

        if(SPI_processed == 0)
            break;

        for(int i = 0; i < SPI_processed; i++)
        {
            HeapTuple tuple = SPI_tuptable->vals[i];

            Datum id = SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);

            lucy_delete(&lucy, DatumGetInt32(id));
        }
    }

    SPI_cursor_close(auditCursor);


    if(unlikely(SPI_exec("delete from " MOLECULES_TABLE " tbl using "
            AUDIT_TABLE " aud where tbl.id = aud.id", 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);


    /*
     * convert new data
     */

    /* TODO: find a better way to obtain bigint[] OID */
    if(unlikely(SPI_exec("select typarray from pg_type where typname = 'int8'", 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    Oid int8arrayOid = DatumGetObjectId(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);

    if(int8arrayOid == 0)
        elog(ERROR, "%s: cannot determine bigint[] oid", __func__);


    SPIPlanPtr moleculesPlan = SPI_prepare("insert into " MOLECULES_TABLE " (id, molecule) values ($1,$2)",
            2, (Oid[]) { INT4OID, BYTEAOID });


    Portal compoundCursor = SPI_cursor_open_with_args(NULL, "select cmp.id, cmp.molfile from " COMPOUNDS_TABLE " cmp, "
            AUDIT_TABLE " aud where cmp.id = aud.id and aud.stored",
            0, NULL, NULL, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);


    VarChar **molfiles = palloc(SYNC_FETCH_SIZE * sizeof(VarChar *));
    LucyLoaderData *data = palloc(SYNC_FETCH_SIZE * sizeof(LucyLoaderData));
    Fingerprint *result = palloc(SYNC_FETCH_SIZE * sizeof(Fingerprint));
    int count = 0;


    while(true)
    {
        SPI_cursor_fetch(compoundCursor, true, SYNC_FETCH_SIZE);

        if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
            elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

        if(SPI_processed == 0)
            break;

        int processed = SPI_processed;
        SPITupleTable *tuptable = SPI_tuptable;

        for(int i = 0; i < processed; i++)
        {
            HeapTuple tuple = tuptable->vals[i];

            Datum molfile = SPI_getbinval(tuple, tuptable->tupdesc, 2, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);

            molfiles[i] = DatumGetVarCharP(molfile);
        }

        java_parse_lucy_data(processed, molfiles, data);

        for(int i = 0; i < processed; i++)
        {
            HeapTuple tuple = SPI_tuptable->vals[i];
            char isNullFlag;

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);
        }

        pthread_t threads[countOfThread];
        ThreadData threadData = {.count = processed, .data = data, .result = result};
        pg_atomic_init_u32(&threadData.possition, 0);

        PG_MEMCONTEXT_BEGIN(memoryContext);
        for(int i = 0; i < countOfThread; i++)
            pthread_create(&threads[i], &attr, lucy_substructure_process_data, (void *) &threadData);
        PG_MEMCONTEXT_END();

        void *status;

        for(int i = 0; i < countOfThread; i++)
            pthread_join(threads[i], &status);

        PG_TRY();
        {
            for(int i = 0; i < processed; i++)
            {
                HeapTuple tuple = tuptable->vals[i];

                if((char *) molfiles[i] != DatumGetPointer(SPI_getbinval(tuple, tuptable->tupdesc, 2, &isNullFlag)))
                    pfree(molfiles[i]);


                Datum id = SPI_getbinval(tuple, tuptable->tupdesc, 1, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);


                if(data[i].error || !fingerprint_is_valid(result[i]))
                {
                    if(data[i].error == NULL)
                        data[i].error = cstring_to_text("fingerprint cannot be generated");

                    char *message = text_to_cstring(data[i].error);
                    elog(NOTICE, "%i: %s", DatumGetInt32(id), message);
                    pfree(message);

                    Datum values[] = {id, PointerGetDatum(data[i].error)};

                    if(SPI_execute_with_args("insert into " MOLECULE_ERRORS_TABLE " (compound, message) values ($1,$2)",
                            2, (Oid[]) { INT4OID, TEXTOID }, values, NULL, false, 0) != SPI_OK_INSERT)
                        elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);

                    pfree(data[i].error);
                }

                if(data[i].molecule == NULL)
                    continue;


                Datum moleculesValues[] = {id, PointerGetDatum(data[i].molecule)};

                if(SPI_execute_plan(moleculesPlan, moleculesValues, NULL, false, 0) != SPI_OK_INSERT)
                    elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

                if(fingerprint_is_valid(result[i]))
                    lucy_add(&lucy, DatumGetInt32(id), result[i]);

                pfree(data[i].molecule);
            }
        }
        PG_CATCH();
        {
            for(int i = 0; i < processed; i++)
                free(result[i].data);

            PG_RE_THROW();
        }
        PG_END_TRY();

        for(int i = 0; i < processed; i++)
            free(result[i].data);

        SPI_freetuptable(tuptable);


        count += processed;

        if(verbose)
            elog(NOTICE, "already processed: %i", count);
    }

    SPI_cursor_close(compoundCursor);

    if(unlikely(SPI_exec("delete from " AUDIT_TABLE, 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);


    lucy_commit(&lucy);
    SPI_finish();
    PG_RETURN_VOID();
}


PG_FUNCTION_INFO_V1(lucy_substructure_optimize_index);
Datum lucy_substructure_optimize_index(PG_FUNCTION_ARGS)
{
    lucy_optimize(&lucy);

    PG_RETURN_BOOL(true);
}
