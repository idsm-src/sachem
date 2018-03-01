#include <postgres.h>
#include <executor/spi.h>
#include <unistd.h>
#include <pthread.h>
#include "molecule.h"
#include "sachem.h"
#include "stats.h"


#define SYNC_FETCH_SIZE           100000
#define MOLECULES_TABLE           "sachem_molecules"
#define ORDER_FILE                "fporder.bin"


typedef struct
{
    pg_atomic_uint32 *possition;
    uint32_t count;
    Stats *stats;
    bytea **molecule;
} ThreadData;


static void *process_data(void *dataPtr)
{
    ThreadData *data = (ThreadData *) dataPtr;
    Molecule molecule;

    while(true)
    {
        bool error = false;

        uint32_t i = pg_atomic_fetch_add_u32(data->possition, 1);

        if(i >= data->count)
            break;

        if(molecule_simple_init(&molecule, VARDATA(data->molecule[i]), malloc))
            if(!stats_add(data->stats, &molecule))
                error = true;

        molecule_simple_free(&molecule, free);

        if(error)
            pthread_exit((void *) 1);
    }

    pthread_exit((void *) 0);
}


PG_FUNCTION_INFO_V1(sachem_generate_fporder);
Datum sachem_generate_fporder(PG_FUNCTION_ARGS)
{
    int32_t limit = PG_GETARG_INT32(0);
    bool verbose = PG_GETARG_BOOL(1);


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    char isNullFlag;


    /*
     * perepare reqired thread data
     */
    int countOfThread = sysconf(_SC_NPROCESSORS_ONLN);

    pthread_attr_t attr;
    pthread_attr_init(&attr);


    ThreadData threadData[countOfThread];

    for(int i = 0; i < countOfThread; i++)
        threadData[i].stats = NULL;

    Stats *stats = NULL;

    PG_TRY();
    {
        if((stats = stats_create()) == NULL)
            elog(ERROR, "%s: stats_create() failed", __func__);

        Portal compoundCursor = SPI_cursor_open_with_args(NULL, "select molecule from " MOLECULES_TABLE,
                0, NULL, NULL, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

        bytea **molecules = palloc(SYNC_FETCH_SIZE * sizeof(bytea *));
        int count = 0;


        while(true)
        {
            SPI_cursor_fetch(compoundCursor, true, SYNC_FETCH_SIZE);

            if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
                elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

            if(SPI_processed == 0)
                break;

            int processed = SPI_processed;
            SPITupleTable *tuptable = SPI_tuptable;

            for(int i = 0; i < processed; i++)
            {
                HeapTuple tuple = tuptable->vals[i];

                Datum mol = SPI_getbinval(tuple, tuptable->tupdesc, 1, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                molecules[i] = DatumGetByteaP(mol);
            }


            pthread_t threads[countOfThread];
            pg_atomic_uint32 possition;
            pg_atomic_init_u32(&possition, 0);

            for(int i = 0; i < countOfThread; i++)
            {
                threadData[i] = (ThreadData){.possition = &possition, .count = processed, .molecule = molecules};

                if((threadData[i].stats = stats_create()) == NULL)
                    elog(ERROR, "%s: stats_create() failed", __func__);
            }

            for(int i = 0; i < countOfThread; i++)
                pthread_create(&threads[i], &attr, process_data, (void *) (threadData + i));

            void *status;

            for(int i = 0; i < countOfThread; i++)
            {
                pthread_join(threads[i], &status);

                if(!stats_merge(stats, threadData[i].stats))
                    elog(ERROR, "%s: stats_merge() failed", __func__);

                stats_delete(threadData[i].stats);
                threadData[i].stats = NULL;
            }

            for(int i = 0; i < processed; i++)
            {
                HeapTuple tuple = tuptable->vals[i];

                if((char *) molecules[i] != DatumGetPointer(SPI_getbinval(tuple, tuptable->tupdesc, 1, &isNullFlag)))
                    pfree(molecules[i]);
            }

            SPI_freetuptable(tuptable);


            count += processed;

            if(verbose)
                elog(NOTICE, "already processed: %i", count);
        }

        SPI_cursor_close(compoundCursor);


        Name database = DatumGetName(DirectFunctionCall1(current_database, 0));
        size_t basePathLength = strlen(DataDir);
        size_t databaseLength = strlen(database->data);

        char *fporderPath = (char *) palloc(basePathLength +  databaseLength + strlen(ORDER_FILE) + 3);
        sprintf(fporderPath, "%s/%s/" ORDER_FILE, DataDir, database->data);

        if(!stats_write(stats, fporderPath, limit))
            elog(ERROR, "%s: stats_write() failed", __func__);

        stats_delete(stats);
    }
    PG_CATCH();
    {
        for(int i = 0; i < countOfThread; i++)
            if(threadData[i].stats != NULL)
                stats_delete(threadData[i].stats);

        if(stats != NULL)
            stats_delete(stats);

        PG_RE_THROW();
    }
    PG_END_TRY();


    SPI_finish();
    PG_RETURN_VOID();
}
