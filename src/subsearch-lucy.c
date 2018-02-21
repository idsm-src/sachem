#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <fmgr.h>
#include <funcapi.h>
#include <utils/array.h>
#include <utils/builtins.h>
#include <utils/memutils.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <postmaster/postmaster.h>
#include <postmaster/fork_process.h>
#include <storage/dsm.h>
#include <storage/pg_shmem.h>
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
#define SYNC_FETCH_SIZE           50000     // per process
#define COMPOUNDS_TABLE           "compounds"
#define MOLECULES_TABLE           "sachem_molecules"
#define MOLECULE_ERRORS_TABLE     "sachem_molecule_errors"
#define AUDIT_TABLE               "sachem_compound_audit"
#define INDEX_TABLE               "sachem_index"


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

    LucyResultSet resultSet;

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


static bool lucyInitialised = false;
static int indexId = -1;
static SPIPlanPtr mainQueryPlan;
static SPIPlanPtr snapshotQueryPlan;
static Lucy lucy;


void lucy_subsearch_init(void)
{
    /* prepare lucy */
    if(unlikely(lucyInitialised == false))
    {
        lucy_init(&lucy);
        lucyInitialised = true;
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

        lucy_set_folder(&lucy, TextDatumGetCString(path));
        indexId = DatumGetInt32(id);
    }

    SPI_finish();
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

        lucy_subsearch_init();

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


                    if(!lucy_is_open(&lucy, &info->resultSet))
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

                        info->resultSet = lucy_search(&lucy, fp, INT32_MAX);


                        if(fp.data != NULL)
                            pfree(fp.data);

                        PG_MEMCONTEXT_END();
                    }


                    int32 *arrayData = (int32 *) ARR_DATA_PTR(info->arrayBuffer);
                    SubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                    int count = lucy_get(&lucy, &info->resultSet, arrayData, FETCH_SIZE);


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
        lucy_fail(&lucy, &info->resultSet);

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


PG_FUNCTION_INFO_V1(lucy_sync_data);
Datum lucy_sync_data(PG_FUNCTION_ARGS)
{
    bool verbose = PG_GETARG_BOOL(0);
    bool optimize = PG_GETARG_BOOL(1);

    MemoryContext memoryContext = CurrentMemoryContext;
    createBasePath();


    /* prepare lucy */
    if(unlikely(lucyInitialised == false))
    {
        lucy_init(&lucy);
        lucyInitialised = true;
    }


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    char isNullFlag;


    /* clone old index */
    int indexNumber = 0;
    char *oldIndexPath = NULL;

    if(unlikely(SPI_exec("select id, path from " INDEX_TABLE, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(SPI_processed != 0)
    {
        if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
            elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

        Datum number = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE) || isNullFlag)
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        indexNumber = DatumGetInt32(number) + 1;


        Datum path = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE) || isNullFlag)
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        oldIndexPath = TextDatumGetCString(path);
    }

    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));
    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);

    int countOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);

    char *indexPath = (char *) palloc(basePathLength +  databaseLength + 64);
    sprintf(indexPath, "%s/%s/lucy-%i", DataDir, database->data, indexNumber);

    char *subindexPath[countOfProcessors];

    for(int p = 0; p < countOfProcessors; p++)
    {
        subindexPath[p] = (char *) palloc(basePathLength +  databaseLength + 64);
        sprintf(subindexPath[p], "%s/%s/lucy-%i.%i", DataDir, database->data, indexNumber, p);
    }


    lucy_link_directory(oldIndexPath, indexPath);
    lucy_set_folder(&lucy, indexPath);

    if(unlikely(SPI_exec("delete from " INDEX_TABLE, 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(SPI_execute_with_args("insert into " INDEX_TABLE " (id, path) values ($1,$2)", 2, (Oid[]) { INT4OID, TEXTOID },
            (Datum[]) {Int32GetDatum(indexNumber), CStringGetTextDatum(indexPath)}, NULL, false, 0) != SPI_OK_INSERT)
        elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);


    lucy_begin(&lucy);

    PG_TRY();
    {
        /*
         * delete unnecessary data
         */

        Portal auditCursor = SPI_cursor_open_with_args(NULL, "select id from " AUDIT_TABLE " where not stored",
                0, NULL, NULL, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

        while(true)
        {
            SPI_cursor_fetch(auditCursor, true, countOfProcessors * SYNC_FETCH_SIZE);

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


        Datum *ids = palloc(countOfProcessors * SYNC_FETCH_SIZE * sizeof(Datum));
        VarChar **molfiles = palloc(countOfProcessors * SYNC_FETCH_SIZE * sizeof(VarChar *));
        LucyLoaderData *data = palloc(countOfProcessors * SYNC_FETCH_SIZE * sizeof(LucyLoaderData));
        int count = 0;


        while(true)
        {
            SPI_cursor_fetch(compoundCursor, true, countOfProcessors * SYNC_FETCH_SIZE);

            if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
                elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

            if(SPI_processed == 0)
                break;

            int processed = SPI_processed;
            SPITupleTable *tuptable = SPI_tuptable;

            for(int i = 0; i < processed; i++)
            {
                HeapTuple tuple = tuptable->vals[i];


                Datum id = SPI_getbinval(tuple, tuptable->tupdesc, 1, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                ids[i] = id;


                Datum molfile = SPI_getbinval(tuple, tuptable->tupdesc, 2, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                molfiles[i] = DatumGetVarCharP(molfile);
            }

            java_parse_lucy_data(processed, molfiles, data);


            int process[countOfProcessors];
            int started = false;

            for(int p = 0; p < countOfProcessors; p++)
            {
                int pid;

                switch((pid = fork_process()))
                {
                    case 0:
                        /* in postmaster child ... */
                        InitPostmasterChild();

                        /* Close the postmaster's sockets */
                        //ClosePostmasterPorts(false);

                        /* Drop our connection to postmaster's shared memory, as well */
                        dsm_detach_all();
                        PGSharedMemoryDetach();

                        PG_TRY();
                        {
                            lucy_set_folder(&lucy, subindexPath[p]);
                            lucy_begin(&lucy);

                            for(int i = p * SYNC_FETCH_SIZE; i < processed && i < (p + 1) * SYNC_FETCH_SIZE; i++)
                            {
                                if(data[i].molecule)
                                {
                                    Molecule molecule;

                                    if(!molecule_simple_init(&molecule, VARDATA(data[i].molecule), malloc))
                                        _exit(1);

                                    Fingerprint result = fingerprint_get(&molecule, &malloc);

                                    if(fingerprint_is_valid(result))
                                        lucy_add(&lucy, DatumGetInt32(ids[i]), result);

                                    molecule_simple_free(&molecule, free);
                                }
                            }

                            lucy_commit(&lucy);
                        }
                        PG_CATCH();
                        {
                            _exit(2);
                        }
                        PG_END_TRY();
                        _exit(0);

                    default:
                        if(pid != -1)
                            started = true;

                        process[p] = pid;
                }
            }

            if(!started)
                elog(ERROR, "index processes failed");


            PG_TRY();
            {
                for(int i = 0; i < processed; i++)
                {
                    HeapTuple tuple = tuptable->vals[i];

                    if((char *) molfiles[i] != DatumGetPointer(SPI_getbinval(tuple, tuptable->tupdesc, 2, &isNullFlag)))
                        pfree(molfiles[i]);


                    if(data[i].error)
                    {
                        if(data[i].error == NULL)
                            data[i].error = cstring_to_text("fingerprint cannot be generated");

                        char *message = text_to_cstring(data[i].error);
                        elog(NOTICE, "%i: %s", DatumGetInt32(ids[i]), message);
                        pfree(message);

                        Datum values[] = {ids[i], PointerGetDatum(data[i].error)};

                        if(SPI_execute_with_args("insert into " MOLECULE_ERRORS_TABLE " (compound, message) values ($1,$2)",
                                2, (Oid[]) { INT4OID, TEXTOID }, values, NULL, false, 0) != SPI_OK_INSERT)
                            elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);

                        pfree(data[i].error);
                    }

                    if(data[i].molecule == NULL)
                        continue;


                    Datum moleculesValues[] = {ids[i], PointerGetDatum(data[i].molecule)};

                    if(SPI_execute_plan(moleculesPlan, moleculesValues, NULL, false, 0) != SPI_OK_INSERT)
                        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

                    pfree(data[i].molecule);
                }
            }
            PG_CATCH();
            {
                int wstatus;

                for(int p = 0; p < countOfProcessors; p++)
                    if(process[p] != -1)
                        waitpid(process[p], &wstatus, 0);

                PG_RE_THROW();
            }
            PG_END_TRY();


            bool error = false;

            for(int p = 0; p < countOfProcessors; p++)
            {
                if(process[p] != -1)
                {
                    int wstatus;

                    waitpid(process[p], &wstatus, 0);

                    if(!WIFEXITED(wstatus) || WEXITSTATUS(wstatus) != 0)
                        error = true;
                }
            }

            if(error)
                elog(ERROR, "an index process failed");


            SPI_freetuptable(tuptable);


            count += processed;

            if(verbose)
                elog(NOTICE, "already processed: %i", count);
        }

        SPI_cursor_close(compoundCursor);

        if(unlikely(SPI_exec("delete from " AUDIT_TABLE, 0) != SPI_OK_DELETE))
            elog(ERROR, "%s: SPI_exec() failed", __func__);


        int process[countOfProcessors];
        int started = false;

        for(int p = 0; p < countOfProcessors; p++)
        {
            int pid;

            switch((pid = fork_process()))
            {
                case 0:
                    /* in postmaster child ... */
                    InitPostmasterChild();

                    /* Close the postmaster's sockets */
                    //ClosePostmasterPorts(false);

                    /* Drop our connection to postmaster's shared memory, as well */
                    dsm_detach_all();
                    PGSharedMemoryDetach();

                    PG_TRY();
                    {
                        lucy_set_folder(&lucy, subindexPath[p]);
                        lucy_begin(&lucy);

                        if(optimize)
                            lucy_optimize(&lucy);

                        lucy_commit(&lucy);
                    }
                    PG_CATCH();
                    {
                        _exit(2);
                    }
                    PG_END_TRY();
                    _exit(0);

                default:
                    if(pid != -1)
                        started = true;

                    process[p] = pid;
            }
        }

        if(!started)
            elog(ERROR, "index processes failed");


        bool error = false;

        for(int p = 0; p < countOfProcessors; p++)
        {
            if(process[p] != -1)
            {
                int wstatus;

                waitpid(process[p], &wstatus, 0);

                if(WIFEXITED(wstatus) && WEXITSTATUS(wstatus) == 0)
                {
                    lucy_add_index(&lucy, subindexPath[p]);
                    lucy_delete_directory(subindexPath[p]);
                }
                else
                {
                    error = true;
                }
            }
        }

        if(error)
            elog(ERROR, "an optimize process failed");


        if(optimize)
            lucy_optimize(&lucy);

        lucy_commit(&lucy);
    }
    PG_CATCH();
    {
        lucy_rollback(&lucy);

        PG_RE_THROW();
    }
    PG_END_TRY();


    SPI_finish();
    PG_RETURN_VOID();
}
