#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <unistd.h>
#include <access/xact.h>
#include <access/parallel.h>
#include <storage/shm_toc.h>
#include <storage/ipc.h>
#include <storage/spin.h>
#include "fporder.h"
#include "molecule.h"
#include "sachem.h"
#include "stats.h"


#define HEADER_KEY                0
#define QUEUE_KEY                 1
#define SYNC_FETCH_SIZE           100000
#define QUEUE_SIZE                1000
#define MOLECULES_TABLE           "sachem_molecules"


typedef struct
{
    slock_t mutex;
    int worker;
    int offset;
    int processed;
    bool verbose;
} WorkerHeader;


void fporder_worker(dsm_segment *seg, shm_toc *toc)
{
    volatile WorkerHeader *header = shm_toc_lookup_key(toc, HEADER_KEY);

    SpinLockAcquire(&header->mutex);
    int worker = header->worker++;
    SpinLockRelease(&header->mutex);

    shm_mq *queue = shm_toc_lookup_key(toc, QUEUE_KEY) + worker * QUEUE_SIZE * sizeof(StatItem);
    shm_mq_set_sender(queue, MyProc);
    shm_mq_handle *out = shm_mq_attach(queue, seg, NULL);


    Stats *stats = stats = stats_create();


    PG_TRY();
    {
        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
            elog(ERROR, "%s: SPI_connect() failed", __func__);

        char isNullFlag;

        SPIPlanPtr queryPlan = SPI_prepare("select molecule from " MOLECULES_TABLE " offset $1 limit $2", 2, (Oid[]) { INT4OID, INT4OID });


        while(true)
        {
            SpinLockAcquire(&header->mutex);
            int offset = header->offset;
            header->offset += SYNC_FETCH_SIZE;
            SpinLockRelease(&header->mutex);


            Datum values[] = { Int32GetDatum(offset), Int32GetDatum(SYNC_FETCH_SIZE)};

            if(unlikely(SPI_execute_plan(queryPlan, values, NULL, true, 0) != SPI_OK_SELECT))
                elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

            if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
                elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

            if(SPI_processed == 0)
                break;

            int processed = SPI_processed;
            SPITupleTable *tuptable = SPI_tuptable;

            for(int i = 0; i < processed; i++)
            {
                CHECK_FOR_INTERRUPTS();

                HeapTuple tuple = tuptable->vals[i];

                Datum mol = SPI_getbinval(tuple, tuptable->tupdesc, 1, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                bytea *data = DatumGetByteaP(mol);

                Molecule molecule;
                molecule_simple_init(&molecule, (uint8_t *) VARDATA(data));

                stats_add(stats, &molecule);

                molecule_simple_free(&molecule);

                if((char *) data != DatumGetPointer(mol))
                    pfree(data);
            }

            SPI_freetuptable(tuptable);

            SpinLockAcquire(&header->mutex);
            header->processed += processed;
            int count = header->processed;
            SpinLockRelease(&header->mutex);

            if(header->verbose)
                elog(NOTICE, "already processed: %i", count);
        }

        SPI_finish();


        StatItem *items;
        size_t count = stats_get_items(stats, &items);

        for(size_t i = 0; i < count; i += QUEUE_SIZE)
        {
            size_t size = count - i;

            if(size > QUEUE_SIZE)
                size = QUEUE_SIZE;

            shm_mq_result result = shm_mq_send(out, size * sizeof(StatItem), items + i, false);

            if (result != SHM_MQ_SUCCESS)
                elog(ERROR, "%s: shm_mq_send() failed", __func__);
        }
    }
    PG_CATCH();
    {
        stats_delete(stats);

        PG_RE_THROW();
    }
    PG_END_TRY();


    shm_mq_detach(queue);
    stats_delete(stats);
}


PG_FUNCTION_INFO_V1(sachem_generate_fporder);
Datum sachem_generate_fporder(PG_FUNCTION_ARGS)
{
    int32_t limit = PG_GETARG_INT32(0);
    bool verbose = PG_GETARG_BOOL(1);

    int countOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);


    Stats *stats = stats_create();

    PG_TRY();
    {
        EnterParallelMode();

        ParallelContext *pcxt = CreateParallelContextForExternalFunction("libsachem", "fporder_worker", countOfProcessors);

        shm_toc_estimate_keys(&pcxt->estimator, 2);
        shm_toc_estimate_chunk(&pcxt->estimator, sizeof(WorkerHeader));
        shm_toc_estimate_chunk(&pcxt->estimator, countOfProcessors * QUEUE_SIZE * sizeof(StatItem));

        InitializeParallelDSM(pcxt);

        WorkerHeader *header = shm_toc_allocate(pcxt->toc, sizeof(WorkerHeader));
        SpinLockInit(&header->mutex);
        header->offset = 0;
        header->processed = 0;
        header->verbose = verbose;
        shm_toc_insert(pcxt->toc, HEADER_KEY, header);

        void *queueBase = shm_toc_allocate(pcxt->toc, countOfProcessors * QUEUE_SIZE * sizeof(StatItem));
        shm_toc_insert(pcxt->toc, QUEUE_KEY, queueBase);

        for(int w = 0; w < countOfProcessors; w++)
            shm_mq_create(queueBase + w * QUEUE_SIZE * sizeof(StatItem), QUEUE_SIZE * sizeof(StatItem));


        LaunchParallelWorkers(pcxt);

        for(int w = 0; w < pcxt->nworkers_launched; w++)
        {
            shm_mq *queue = queueBase + w * QUEUE_SIZE * sizeof(StatItem);
            shm_mq_set_receiver(queue, MyProc);
            shm_mq_handle *in = shm_mq_attach(queue, pcxt->seg, NULL);

            while(true)
            {
                Size bytes;
                StatItem *items;
                shm_mq_result result = shm_mq_receive(in, &bytes, (void *) &items, false);

                if(result == SHM_MQ_SUCCESS)
                    stats_merge(stats, items, bytes / sizeof(StatItem));
                else if(result != SHM_MQ_DETACHED)
                    elog(ERROR, "%s: shm_mq_receive() failed", __func__);
                else
                    break;
            }

            shm_mq_detach(queue);
        }

        WaitForParallelWorkersToFinish(pcxt);
        DestroyParallelContext(pcxt);
        ExitParallelMode();


        char *fporderPath = get_file_path(ORDER_FILE);
        stats_write(stats, fporderPath, limit);
    }
    PG_CATCH();
    {
        stats_delete(stats);

        PG_RE_THROW();
    }
    PG_END_TRY();

    stats_delete(stats);
    PG_RETURN_VOID();
}
