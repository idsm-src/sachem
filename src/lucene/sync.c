#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <access/xact.h>
#include <access/parallel.h>
#include <storage/shm_toc.h>
#include <storage/spin.h>
#include <unistd.h>
#include "common.h"
#include "molecule.h"
#include "sachem.h"
#include "lucene.h"
#include "java/parse.h"
#include "fingerprints/fingerprint.h"


#define SYNC_FETCH_SIZE         100000
#define SUBOPTIMIZE_PROCESSES   4
#define HEADER_KEY              0
#define ID_TABLE_KEY            1
#define IDEX_KEY_OFFSET         2
#define MOLECULE_KEY_OFFSET     10000


typedef struct
{
    slock_t mutex;
    int attachedWorkers;
    int moleculeCount;
    int moleculePosition;
} IndexWorkerHeader;


typedef struct
{
    slock_t mutex;
    int subindexCount;
    int subindexPosition;
} OptimizeWorkerHeader;


static bool javaInitialized = false;
static bool luceneInitialised = false;
static Lucene lucene;


void lucene_index_worker(dsm_segment *seg, shm_toc *toc)
{
    volatile IndexWorkerHeader *header = shm_toc_lookup_key(toc, HEADER_KEY);
    SpinLockAcquire(&header->mutex);
    int workerNumber = header->attachedWorkers++;
    SpinLockRelease(&header->mutex);

    int *ids = shm_toc_lookup_key(toc, ID_TABLE_KEY);
    char *indexPath = shm_toc_lookup_key(toc, IDEX_KEY_OFFSET + workerNumber);

    lucene_init(&lucene);
    lucene_set_folder(&lucene, indexPath);
    lucene_begin(&lucene);

    PG_TRY();
    {
        while(true)
        {
            SpinLockAcquire(&header->mutex);
            int position = header->moleculePosition++;
            SpinLockRelease(&header->mutex);

            if(position >= header->moleculeCount)
                break;

            uint8_t *data = shm_toc_lookup_key(toc, MOLECULE_KEY_OFFSET + position);

            Molecule molecule;

            molecule_simple_init(&molecule, data);

            IntegerFingerprint result = integer_fingerprint_get(&molecule);
            lucene_add(&lucene, ids[position], result);

            integer_fingerprint_free(result);
            molecule_simple_free(&molecule);
        }

        lucene_commit(&lucene);
    }
    PG_CATCH();
    {
        lucene_rollback(&lucene);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


void lucene_optimize_worker(dsm_segment *seg, shm_toc *toc)
{
    volatile OptimizeWorkerHeader *header = shm_toc_lookup_key(toc, HEADER_KEY);

    lucene_init(&lucene);

    while(true)
    {
        SpinLockAcquire(&header->mutex);
        int position = header->subindexPosition++;
        SpinLockRelease(&header->mutex);

        if(position >= header->subindexCount)
            break;

        char *indexPath = shm_toc_lookup_key(toc, IDEX_KEY_OFFSET + position);

        lucene_set_folder(&lucene, indexPath);
        lucene_begin(&lucene);

        PG_TRY();
        {
            lucene_optimize(&lucene);
            lucene_commit(&lucene);
        }
        PG_CATCH();
        {
            lucene_rollback(&lucene);

            PG_RE_THROW();
        }
        PG_END_TRY();
    }
}


PG_FUNCTION_INFO_V1(lucene_sync_data);
Datum lucene_sync_data(PG_FUNCTION_ARGS)
{
    if(unlikely(javaInitialized == false))
    {
        java_parse_init();
        javaInitialized = true;
    }


    bool verbose = PG_GETARG_BOOL(0);
    bool optimize = PG_GETARG_BOOL(1);

    createBasePath();


    /* prepare lucene */
    if(unlikely(luceneInitialised == false))
    {
        lucene_init(&lucene);
        luceneInitialised = true;
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
    sprintf(indexPath, "%s/%s/lucene-%i", DataDir, database->data, indexNumber);

    char *subindexPath[countOfProcessors];

    for(int p = 0; p < countOfProcessors; p++)
    {
        subindexPath[p] = (char *) palloc(basePathLength +  databaseLength + 64);
        sprintf(subindexPath[p], "%s/%s/lucene-%i.%i", DataDir, database->data, indexNumber, p);
    }


    lucene_link_directory(oldIndexPath, indexPath);
    lucene_set_folder(&lucene, indexPath);

    if(unlikely(SPI_exec("delete from " INDEX_TABLE, 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(SPI_execute_with_args("insert into " INDEX_TABLE " (id, path) values ($1,$2)", 2, (Oid[]) { INT4OID, TEXTOID },
            (Datum[]) {Int32GetDatum(indexNumber), CStringGetTextDatum(indexPath)}, NULL, false, 0) != SPI_OK_INSERT)
        elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);


    lucene_begin(&lucene);
    int subindexCount = 0;


    PG_TRY();
    {
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

                lucene_delete(&lucene, DatumGetInt32(id));
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


        Datum *ids = palloc(SYNC_FETCH_SIZE * sizeof(Datum));
        VarChar **molfiles = palloc(SYNC_FETCH_SIZE * sizeof(VarChar *));
        LoaderData *data = palloc(SYNC_FETCH_SIZE * sizeof(LoaderData));
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


                Datum id = SPI_getbinval(tuple, tuptable->tupdesc, 1, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                ids[i] = id;


                Datum molfile = SPI_getbinval(tuple, tuptable->tupdesc, 2, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                molfiles[i] = DatumGetVarCharP(molfile);
            }

            java_parse_data(processed, molfiles, data);


            EnterParallelMode();

            ParallelContext *pcxt = CreateParallelContextForExternalFunction("libsachem", "lucene_index_worker", countOfProcessors);

            shm_toc_estimate_keys(&pcxt->estimator, 2 + countOfProcessors + processed);
            shm_toc_estimate_chunk(&pcxt->estimator, sizeof(IndexWorkerHeader));
            shm_toc_estimate_chunk(&pcxt->estimator, processed * sizeof(int));

            for(int i = 0; i < countOfProcessors; i++)
                shm_toc_estimate_chunk(&pcxt->estimator, strlen(subindexPath[i]) + 1);

            for(int i = 0; i < processed; i++)
                if(data[i].molecule != NULL)
                    shm_toc_estimate_chunk(&pcxt->estimator, VARSIZE(data[i].molecule) - VARHDRSZ);

            InitializeParallelDSM(pcxt);

            int *idValues = shm_toc_allocate(pcxt->toc, processed * sizeof(int));
            shm_toc_insert(pcxt->toc, ID_TABLE_KEY, idValues);

            for(int i = 0; i < countOfProcessors; i++)
            {
                void *path = shm_toc_allocate(pcxt->toc, strlen(subindexPath[i]) + 1);
                strcpy(path, subindexPath[i]);
                shm_toc_insert(pcxt->toc, IDEX_KEY_OFFSET + i, path);
            }


            int valid = 0;

            for(int i = 0; i < processed; i++)
            {
                if(likely(data[i].molecule != NULL))
                {
                    void *molecule = shm_toc_allocate(pcxt->toc, VARSIZE(data[i].molecule) - VARHDRSZ);
                    memcpy(molecule, VARDATA(data[i].molecule), VARSIZE(data[i].molecule) - VARHDRSZ);
                    shm_toc_insert(pcxt->toc, MOLECULE_KEY_OFFSET + valid, molecule);

                    idValues[valid] = DatumGetInt32(ids[i]);
                    valid++;
                }
            }

            IndexWorkerHeader *header = shm_toc_allocate(pcxt->toc, sizeof(IndexWorkerHeader));
            SpinLockInit(&header->mutex);
            header->attachedWorkers = 0;
            header->moleculePosition = 0;
            header->moleculeCount = valid;
            shm_toc_insert(pcxt->toc, HEADER_KEY, header);

            LaunchParallelWorkers(pcxt);

            if(subindexCount < pcxt->nworkers_launched)
                subindexCount = pcxt->nworkers_launched;

            WaitForParallelWorkersToFinish(pcxt);
            DestroyParallelContext(pcxt);
            ExitParallelMode();


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


            SPI_freetuptable(tuptable);

            count += processed;

            if(verbose)
                elog(NOTICE, "already processed: %i", count);
        }

        SPI_cursor_close(compoundCursor);

        if(unlikely(SPI_exec("delete from " AUDIT_TABLE, 0) != SPI_OK_DELETE))
            elog(ERROR, "%s: SPI_exec() failed", __func__);


        if(optimize)
        {
            EnterParallelMode();

            ParallelContext *pcxt = CreateParallelContextForExternalFunction("libsachem", "lucene_optimize_worker",
                    countOfProcessors > SUBOPTIMIZE_PROCESSES ? SUBOPTIMIZE_PROCESSES : countOfProcessors);

            shm_toc_estimate_keys(&pcxt->estimator, 1 + countOfProcessors);
            shm_toc_estimate_chunk(&pcxt->estimator, sizeof(OptimizeWorkerHeader));

            for(int i = 0; i < countOfProcessors; i++)
                shm_toc_estimate_chunk(&pcxt->estimator, strlen(subindexPath[i]) + 1);

            InitializeParallelDSM(pcxt);

            for(int i = 0; i < countOfProcessors; i++)
            {
                void *path = shm_toc_allocate(pcxt->toc, strlen(subindexPath[i]) + 1);
                strcpy(path, subindexPath[i]);
                shm_toc_insert(pcxt->toc, IDEX_KEY_OFFSET + i, path);
            }

            OptimizeWorkerHeader *header = shm_toc_allocate(pcxt->toc, sizeof(OptimizeWorkerHeader));
            SpinLockInit(&header->mutex);
            header->subindexPosition = 0;
            header->subindexCount = subindexCount;
            shm_toc_insert(pcxt->toc, HEADER_KEY, header);

            LaunchParallelWorkers(pcxt);
            WaitForParallelWorkersToFinish(pcxt);
            DestroyParallelContext(pcxt);
            ExitParallelMode();
        }

        for(int p = 0; p < subindexCount; p++)
        {
            lucene_add_index(&lucene, subindexPath[p]);
            lucene_delete_directory(subindexPath[p]);
        }

        if(optimize)
            lucene_optimize(&lucene);

        lucene_commit(&lucene);
    }
    PG_CATCH();
    {
        lucene_rollback(&lucene);

        PG_RE_THROW();
    }
    PG_END_TRY();


    SPI_finish();
    PG_RETURN_VOID();
}
