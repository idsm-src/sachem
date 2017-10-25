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
#include "pgchem.h"
#include "subsearch.h"
#include "fpsearch/fpsearch-lucy.h"


#define SHOW_STATS                0
#define FETCH_SIZE                1000
#define MOLECULES_TABLE           "orchem_molecules"


typedef struct
{
    int32_t topN;
    bool strictStereo;
    bool exact;

    BitSet resultMask;
    int32_t foundResults;

    SubstructureQueryData *queryData;
    int queryDataCount;
    int queryDataPosition;

    void *lucySearch;

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
static SPIPlanPtr mainQueryPlan;
static void *fplucy;
static pthread_mutex_t indexMutex;


void subsearch_lucy_module_init(void)
{
    /* prepare lucy */
    fplucy_initialize(&fplucy, getFilePath("lucy"), DATADIR "/fporder.txt");

    if(unlikely(fplucy == NULL))
        elog(ERROR, "subsearch-lucy module: lucy initialization failed");


    /* prepare query plan */
    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "subsearch-lucy module: SPI_connect() failed");

    mainQueryPlan = SPI_prepare("select id, seqid, atoms, bonds from " MOLECULES_TABLE " where seqid = any($1)", 1, (Oid[]) { INT4ARRAYOID });

    if(unlikely(mainQueryPlan == NULL))
        elog(ERROR, "subsearch-lucy module: SPI_prepare_cursor() failed");

    if(unlikely(SPI_keepplan(mainQueryPlan) == SPI_ERROR_ARGUMENT))
        elog(ERROR, "subsearch-lucy module: SPI_keepplan() failed");


    SPI_finish();
    initialised = true;
}


void subsearch_lucy_module_finish(void)
{
    if(fplucy != NULL)
        fplucy_close(fplucy);

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
            elog(ERROR, "subsearch-lucy module is not properly initialized");

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        text *type = PG_GETARG_TEXT_P(1);
        int32_t topN = PG_GETARG_INT32(2);
        bool strictStereo = PG_GETARG_BOOL(3);
        bool exact = PG_GETARG_BOOL(4);
        bool tautomers = PG_GETARG_BOOL(5);
        char *typeStr = text_to_cstring(type);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();


        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
             elog(ERROR, "subsearch-lucy module: SPI_connect() failed");

        connected = true;

        if(unlikely(SPI_execute("select max(seqid) + 1 from " MOLECULES_TABLE, true, FETCH_ALL) != SPI_OK_SELECT))
            elog(ERROR, "subsearch-lucy module: SPI_execute() failed");

        if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
            elog(ERROR, "subsearch-lucy module: SPI_execute() failed");

        char isNullFlag;
        int64_t moleculeCount = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "subsearch-lucy module: SPI_getbinval() failed");

        SPI_freetuptable(SPI_tuptable);


        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);

        SubstructureSearchData *info = (SubstructureSearchData *) palloc(sizeof(SubstructureSearchData));
        funcctx->user_fctx = info;

        info->topN = topN;
        info->strictStereo = strictStereo;
        info->exact = exact;

        info->queryDataCount = java_parse_substructure_query(&info->queryData, VARDATA(query), VARSIZE(query) - VARHDRSZ, typeStr, tautomers);

        PG_FREE_IF_COPY(query, 0);
        PG_FREE_IF_COPY(type, 1);
        pfree(typeStr);

        info->queryDataPosition = -1;
        info->lucySearch = NULL;
        info->table = NULL;
        info->tableRowCount = -1;
        info->tableRowPosition = -1;
        info->foundResults = 0;

        bitset_init_empty(&info->resultMask, moleculeCount);

        info->isomorphismContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx, "subsearch-lucy isomorphism context",
                ALLOCSET_DEFAULT_MINSIZE, ALLOCSET_DEFAULT_INITSIZE, ALLOCSET_DEFAULT_MAXSIZE);
        info->targetContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx, "subsearch-lucy target context",
                ALLOCSET_DEFAULT_MINSIZE, ALLOCSET_DEFAULT_INITSIZE, ALLOCSET_DEFAULT_MAXSIZE);

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
                    molecule_init(&info->queryMolecule, data->atomLength, data->atoms, data->bondLength, data->bonds, data->restH, !info->exact);
                    vf2state_init(&info->vf2state, &info->queryMolecule, info->strictStereo, info->exact);
                    PG_MEMCONTEXT_END();

                    info->lucySearch = fplucy_search(fplucy, &info->queryMolecule, INT32_MAX);
                }


                int32 *arrayData = (int32 *) ARR_DATA_PTR(info->arrayBuffer);
                SubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                int count = fplucy_search_fillbuf(fplucy, info->lucySearch, arrayData, FETCH_SIZE);


                if(unlikely(count == 0))
                {
                    fplucy_search_finish(fplucy, info->lucySearch);
                    info->lucySearch = NULL;

                    continue;
                }


                *(ARR_DIMS(info->arrayBuffer)) = count;
                SET_VARSIZE(info->arrayBuffer, count * sizeof(int32) + ARR_OVERHEAD_NONULLS(1));

                Datum values[] = { PointerGetDatum(info->arrayBuffer)};


                if(unlikely(!connected && SPI_connect() != SPI_OK_CONNECT))
                     elog(ERROR, "subsearch-lucy module: SPI_connect() failed");

                connected = true;


                if(unlikely(SPI_execute_plan(mainQueryPlan, values, NULL, true, 0) != SPI_OK_SELECT))
                    elog(ERROR, "subsearch-lucy module: SPI_execute_plan() failed");

                if(unlikely(SPI_processed != count || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 4))
                    elog(ERROR, "subsearch-lucy module: SPI_execute_plan() failed");

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
                elog(ERROR, "subsearch-lucy module: SPI_getbinval() failed");


            int32 seqid = DatumGetInt32(SPI_getbinval(tuple, tupdesc, 2, &isNullFlag));

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch-lucy module: SPI_getbinval() failed");

            if(bitset_get(&info->resultMask, seqid))
                continue;


            Datum atoms = SPI_getbinval(tuple, tupdesc, 3, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch-lucy module: SPI_getbinval() failed");


            Datum bonds = SPI_getbinval(tuple, tupdesc, 4, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch-lucy module: SPI_getbinval() failed");

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

            bool match;

            PG_MEMCONTEXT_BEGIN(info->targetContext);
            Molecule target;
            molecule_init(&target, VARSIZE(atomsData) - VARHDRSZ, VARDATA(atomsData), VARSIZE(bondsData) - VARHDRSZ, VARDATA(bondsData), NULL, false);
            match = vf2state_match(&info->vf2state, &target);
            PG_MEMCONTEXT_END();
            MemoryContextReset(info->targetContext);

            if(match)
            {
                bitset_set(&info->resultMask, seqid);
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


void *lucy_substructure_process_spi_table(void* idx)
{
    pthread_mutex_lock(&indexMutex);
    MemoryContext indexContext = AllocSetContextCreate(CurrentMemoryContext, "index thread context",
            ALLOCSET_DEFAULT_MINSIZE, ALLOCSET_DEFAULT_INITSIZE, ALLOCSET_DEFAULT_MAXSIZE);
    pthread_mutex_unlock(&indexMutex);


    while(true)
    {
        pthread_mutex_lock(&indexMutex);
        MemoryContextReset(indexContext);

        volatile int *index = (volatile int*) idx;
        int i = *index;

        if(i >= SPI_processed)
            break;

        HeapTuple tuple = SPI_tuptable->vals[i];
        char isNullFlag;


        Datum seqid = SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "lucy_substructure_create_index: SPI_getbinval() failed");


        Datum atoms = SPI_getbinval(tuple, SPI_tuptable->tupdesc, 2, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "lucy_substructure_create_index: SPI_getbinval() failed");


        Datum bonds = SPI_getbinval(tuple, SPI_tuptable->tupdesc, 3, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "lucy_substructure_create_index: SPI_getbinval() failed");


        Molecule molecule;

        PG_MEMCONTEXT_BEGIN(indexContext);
        bytea *atomsData = DatumGetByteaP(atoms);
        bytea *bondsData = DatumGetByteaP(bonds);
        molecule_init(&molecule, VARSIZE(atomsData) - VARHDRSZ, VARDATA(atomsData), VARSIZE(bondsData) - VARHDRSZ, VARDATA(bondsData), NULL, false);
        PG_MEMCONTEXT_END();

        (*index)++;
        pthread_mutex_unlock(&indexMutex);

        fplucy_add_mol(fplucy, DatumGetInt32(seqid), &molecule);
    }


    pthread_mutex_unlock(&indexMutex);
    pthread_exit((void*) 0);
}


PG_FUNCTION_INFO_V1(lucy_substructure_create_index);
Datum lucy_substructure_create_index(PG_FUNCTION_ARGS)
{
    MemoryContext indexContext = AllocSetContextCreate(CurrentMemoryContext, "lucy indexing context",
            ALLOCSET_DEFAULT_MINSIZE, ALLOCSET_DEFAULT_INITSIZE, ALLOCSET_DEFAULT_MAXSIZE);


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
         elog(ERROR, "lucy_substructure_create_index: SPI_connect() failed");


    Portal cursor = SPI_cursor_open_with_args(NULL, "select seqid, atoms, bonds from " MOLECULES_TABLE,
            0, NULL, NULL, NULL, true, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

    if(unlikely(cursor == NULL))
            elog(ERROR, "lucy_substructure_create_index: SPI_cursor_open_with_args() failed");


    int countOfThread = sysconf(_SC_NPROCESSORS_ONLN);
    int processed = 0;

    while(true)
    {
        SPI_cursor_fetch(cursor, true, 10000);

        if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 3))
            elog(ERROR, "lucy_substructure_create_index: SPI_cursor_fetch() failed");

        if(SPI_processed == 0)
            break;


        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

        int idx = 0;
        pthread_t threads[countOfThread];
        pthread_mutex_init(&indexMutex, NULL);

        for(int i = 0; i < countOfThread; i++)
            pthread_create(&threads[i], &attr, lucy_substructure_process_spi_table, (void *) &idx);

        pthread_attr_destroy(&attr);

        void *status;

        for(int i = 0; i < countOfThread; i++)
            pthread_join(threads[i], &status);

        pthread_mutex_destroy(&indexMutex);

        processed += SPI_processed;
        elog(NOTICE, "already processed: %i", processed);

        SPI_freetuptable(SPI_tuptable);
    }

    SPI_finish();

    PG_RETURN_BOOL(true);
}
