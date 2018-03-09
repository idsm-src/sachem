#include <postgres.h>
#include <executor/spi.h>
#include <stdbool.h>
#include <unistd.h>
#include "molindex.h"
#include "sachem.h"


#define FETCH_SIZE              100000
#define MOLECULES_TABLE         "sachem_molecules"


void sachem_generate_molecule_index(int indexNumber, bool useSeqId)
{
    char *indexFilePath = get_index_path(MOLECULE_INDEX_PREFIX, MOLECULE_INDEX_SUFFIX, indexNumber);

    int indexFd = open(indexFilePath, O_EXCL | O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);

    if(indexFd == -1)
        elog(ERROR, "%s: open() failed", __func__);


    PG_TRY();
    {
        if(unlikely(SPI_execute("select max(id) + 1 from " MOLECULES_TABLE, false, FETCH_ALL) != SPI_OK_SELECT))
            elog(ERROR, "%s: SPI_execute() failed", __func__);

        if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
            elog(ERROR, "%s: SPI_execute() failed", __func__);

        char isNullFlag;
        uint64_t moleculeCountDatum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        uint64_t moleculeCount = DatumGetInt64(moleculeCountDatum);

        SPI_freetuptable(SPI_tuptable);


        if(lseek(indexFd, (moleculeCount + 1) * sizeof(uint64_t), SEEK_SET) == -1)
            elog(ERROR, "%s: lseek() failed", __func__);


        uint64_t *offsetTable = palloc_extended(moleculeCount * sizeof(uint64_t), MCXT_ALLOC_HUGE);

        for(int i = 0; i < moleculeCount; i++)
            offsetTable[i] = (uint64_t) -1;


        char *query = useSeqId ? "select seqid, molecule, id from " MOLECULES_TABLE " order by seqid" :
                "select id, molecule from " MOLECULES_TABLE " order by id";

        Portal moleculeCursor = SPI_cursor_open_with_args(NULL, query, 0, NULL, NULL, NULL, false,
                CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);


        uint64_t offset = 0;

        while(true)
        {
            SPI_cursor_fetch(moleculeCursor, true, FETCH_SIZE);

            if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2 + useSeqId))
                elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

            if(SPI_processed == 0)
                break;


            for(size_t i = 0; i < SPI_processed; i++)
            {
                HeapTuple tuple = SPI_tuptable->vals[i];

                int id = DatumGetInt32(SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag));

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                Datum moleculeDatum = SPI_getbinval(tuple, SPI_tuptable->tupdesc, 2, &isNullFlag);

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                offsetTable[id] = offset;


                if(useSeqId)
                {
                    int32_t realId = DatumGetInt32(SPI_getbinval(tuple, SPI_tuptable->tupdesc, 3, &isNullFlag));

                    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                        elog(ERROR, "%s: SPI_getbinval() failed", __func__);

                    if(write(indexFd, &realId, sizeof(int32_t)) != sizeof(int32_t))
                        elog(ERROR, "%s: write() failed", __func__);

                    offset += sizeof(int32_t);
                }


                bytea *moleculeData = DatumGetByteaP(moleculeDatum);
                void *data = VARDATA(moleculeData);
                size_t size = VARSIZE(moleculeData) - VARHDRSZ;

                if(write(indexFd, data, size) != size)
                    elog(ERROR, "%s: write() failed", __func__);

                if((void *) moleculeData != DatumGetPointer(moleculeDatum))
                    pfree(moleculeData);

                offset += size;


                if(useSeqId)
                {
                    size_t padding = sizeof(int32_t) - size % sizeof(int32_t);

                    if(padding < sizeof(int32_t))
                    {
                        int32_t zero = 0;

                        if(write(indexFd, &zero, padding) != padding)
                            elog(ERROR, "%s: write() failed", __func__);

                        offset += padding;
                    }
                }
            }

            SPI_freetuptable(SPI_tuptable);
        }

        SPI_cursor_close(moleculeCursor);


        if(lseek(indexFd, 0, SEEK_SET) == -1)
            elog(ERROR, "%s: lseek() failed", __func__);

        if(write(indexFd, &moleculeCount, sizeof(uint64_t)) != sizeof(uint64_t))
            elog(ERROR, "%s: write() failed", __func__);


        int32_t start = 0;

        while(true)
        {
            int32_t skip = 0;

            while(start + skip < moleculeCount && offsetTable[start + skip] == (uint64_t) -1)
                skip++;

            if(lseek(indexFd, skip * sizeof(uint64_t), SEEK_CUR) == -1)
                elog(ERROR, "%s: lseek() failed", __func__);

            start += skip;


            int32_t size = 0;

            while(start + size < moleculeCount && offsetTable[start + size] != (uint64_t) -1)
                size++;

            if(size == 0)
                break;

            if(write(indexFd, offsetTable + start, size * sizeof(uint64_t)) != size * sizeof(uint64_t))
                elog(ERROR, "%s: write() failed", __func__);

            start += size;
        }


        if(close(indexFd) != 0)
            elog(ERROR, "%s: close() failed", __func__);
    }
    PG_CATCH();
    {
        unlink(indexFilePath);

        PG_RE_THROW();
    }
    PG_END_TRY();
}
