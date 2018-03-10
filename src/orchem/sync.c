#include <postgres.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <stdbool.h>
#include <unistd.h>
#include "bitset.h"
#include "common.h"
#include "molindex.h"
#include "sachem.h"
#include "java/orchem.h"


#define SYNC_FETCH_SIZE         100000


static bool javaInitialized = false;


PG_FUNCTION_INFO_V1(orchem_sync_data);
Datum orchem_sync_data(PG_FUNCTION_ARGS)
{
    if(unlikely(javaInitialized == false))
    {
        java_orchem_init();
        javaInitialized = true;
    }


    bool verbose = PG_GETARG_BOOL(0);


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    char isNullFlag;


    /*
     * delete unnecessary data
     */

    if(unlikely(SPI_exec("delete from " MOLECULES_TABLE " tbl using "
            AUDIT_TABLE " aud where tbl.id = aud.id", 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(unlikely(SPI_exec("delete from " MOLECULE_COUNTS_TABLE " tbl using "
            AUDIT_TABLE " aud where tbl.id = aud.id", 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(unlikely(SPI_exec("delete from " FINGERPRINT_TABLE " tbl using "
            AUDIT_TABLE " aud where tbl.id = aud.id", 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);


    /*
     * compoute count of index items
     */

    if(unlikely(SPI_execute("select coalesce(max(seqid)+1, 0), count(*) from "
            MOLECULES_TABLE, false, FETCH_ALL) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute() failed", __func__);

    if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2)
        elog(ERROR, "%s: SPI_execute() failed", __func__);

    int32_t seqidCount = DatumGetInt32(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);

    int32_t moleculesCount = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);


    if(unlikely(SPI_execute("select count(*) from " AUDIT_TABLE " where stored", false, FETCH_ALL) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute() failed", __func__);

    if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
        elog(ERROR, "%s: SPI_execute() failed", __func__);

    int32_t auditCount = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);


    int32_t maxIndexSize = seqidCount > moleculesCount + auditCount ? seqidCount : moleculesCount + auditCount;


    /*
     * get available seqid bitmap
     */

    Portal seqidCursor = SPI_cursor_open_with_args(NULL, "select seqid from " MOLECULES_TABLE,
            0, NULL, NULL, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

    BitSet seqidSet;
    bitset_init_setted(&seqidSet, maxIndexSize);

    while(true)
    {
        SPI_cursor_fetch(seqidCursor, true, 100000);

        if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
            elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

        if(SPI_processed == 0)
            break;


        for(size_t i = 0; i < SPI_processed; i++)
        {
            int32_t seqid = DatumGetInt32(SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &isNullFlag));

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);

            bitset_unset(&seqidSet, seqid);
        }
    }

    SPI_cursor_close(seqidCursor);


    /*
     * load original fingerprint data
     */

    BitSet bitmap[FP_SIZE];

    for(int i = 0; i < FP_SIZE; i++)
        bitset_init_empty(bitmap + i, maxIndexSize);


    Portal fpCursor = SPI_cursor_open_with_args(NULL, "select mt.seqid, ft.fp from " FINGERPRINT_TABLE " ft, "
            MOLECULES_TABLE " mt where mt.id = ft.id", 0, NULL, NULL, NULL, false,
            CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

    while(true)
    {
        SPI_cursor_fetch(fpCursor, true, 100000);

        if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
            elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

        if(SPI_processed == 0)
            break;


        for(size_t i = 0; i < SPI_processed; i++)
        {
            HeapTuple tuple = SPI_tuptable->vals[i];

            int seqid = DatumGetInt32(SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag));

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);

            Datum fpDatum = SPI_getbinval(tuple, SPI_tuptable->tupdesc, 2, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);

            ArrayType *fpArray = DatumGetArrayTypeP(fpDatum);

            if(ARR_NDIM(fpArray))
            {
                int length = ARR_DIMS(fpArray)[0];
                uint64_t *data = (uint64_t *) ARR_DATA_PTR(fpArray);

                BitSet fp;
                bitset_init(&fp, data, length);

                for(int j = bitset_next_set_bit(&fp, 1); j >= 0 && j <= FP_SIZE; j = bitset_next_set_bit(&fp, j + 1))
                    bitset_set(bitmap + j - 1, seqid);
            }
        }

        SPI_freetuptable(SPI_tuptable);
    }

    SPI_cursor_close(fpCursor);


#if USE_COUNT_FINGERPRINT
    /*
     * load original count data
     */

    int16 (*counts)[COUNTS_SIZE] = palloc_extended(maxIndexSize * COUNTS_SIZE * sizeof(int16), MCXT_ALLOC_HUGE | MCXT_ALLOC_ZERO);

    Portal countCursor = SPI_cursor_open_with_args(NULL, "select mt.seqid, ct.counts from " MOLECULE_COUNTS_TABLE " ct, "
            MOLECULES_TABLE " mt where ct.id = mt.id", 0, NULL, NULL, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);


    while(true)
    {
        SPI_cursor_fetch(countCursor, true, 100000);

        if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
            elog(ERROR, "%s: SPI_cursor_fetch() failed", __func__);

        if(SPI_processed == 0)
            break;


        for(size_t i = 0; i < SPI_processed; i++)
        {
            HeapTuple tuple = SPI_tuptable->vals[i];

            int seqid = DatumGetInt32(SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag));

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);

            Datum countsDatum = SPI_getbinval(tuple, SPI_tuptable->tupdesc, 2, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);

            ArrayType *countsArray = DatumGetArrayTypeP(countsDatum);

            if(ARR_NDIM(countsArray) != 1 && ARR_DIMS(countsArray)[0] != COUNTS_SIZE)
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);

            for(int j = 0; j < COUNTS_SIZE; j++)
                counts[seqid][j] = ((int16 *) ARR_DATA_PTR(countsArray))[j];
        }

        SPI_freetuptable(SPI_tuptable);
    }

    SPI_cursor_close(countCursor);
#endif


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


    SPIPlanPtr moleculesPlan = SPI_prepare("insert into " MOLECULES_TABLE " (seqid, id, molecule) values ($1,$2,$3)",
            3, (Oid[]) { INT4OID, INT4OID, BYTEAOID });

    SPIPlanPtr countsPlan = SPI_prepare("insert into " MOLECULE_COUNTS_TABLE " (id, counts) values ($1,$2)",
            2, (Oid[]) { INT4OID, INT2ARRAYOID });

    SPIPlanPtr fingerprintPlan = SPI_prepare("insert into " FINGERPRINT_TABLE " (id, bit_count, fp) values ($1,$2,$3)",
            3, (Oid[]) { INT4OID, INT2OID, int8arrayOid });


    Portal compoundCursor = SPI_cursor_open_with_args(NULL, "select cmp.id, cmp.molfile from " COMPOUNDS_TABLE " cmp, "
            AUDIT_TABLE " aud where cmp.id = aud.id and aud.stored",
            0, NULL, NULL, NULL, false, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);


    VarChar **molfiles = palloc(SYNC_FETCH_SIZE * sizeof(VarChar *));
    OrchemLoaderData *data = palloc(SYNC_FETCH_SIZE * sizeof(OrchemLoaderData));
    int currentSeqid = -1;
    int notIndexed = 0;
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


        java_orchem_parse_data(processed, molfiles, data);


        for(int i = 0; i < processed; i++)
        {
            HeapTuple tuple = tuptable->vals[i];

            if((char *) molfiles[i] != DatumGetPointer(SPI_getbinval(tuple, tuptable->tupdesc, 2, &isNullFlag)))
                pfree(molfiles[i]);


            Datum id = SPI_getbinval(tuple, tuptable->tupdesc, 1, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);


            if(data[i].error)
            {
                char *message = text_to_cstring(data[i].error);
                elog(NOTICE, "%i: %s", DatumGetInt32(id), message);
                pfree(message);

                Datum values[] = {id, PointerGetDatum(data[i].error)};

                if(SPI_execute_with_args("insert into " MOLECULE_ERRORS_TABLE " (compound, message) values ($1,$2)",
                        2, (Oid[]) { INT4OID, TEXTOID }, values, NULL, false, 0) != SPI_OK_INSERT)
                    elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);

                pfree(data[i].error);
            }

            if(data[i].molecule == NULL || data[i].counts == NULL || data[i].fp == NULL)
            {
                notIndexed++;
                continue;
            }


            currentSeqid = bitset_next_set_bit(&seqidSet, currentSeqid + 1);
            Datum moleculesValues[] = {Int32GetDatum(currentSeqid), id, PointerGetDatum(data[i].molecule)};

            if(SPI_execute_plan(moleculesPlan, moleculesValues, NULL, false, 0) != SPI_OK_INSERT)
                elog(ERROR, "%s: SPI_execute_plan() failed", __func__);


#if USE_COUNT_FINGERPRINT
            Datum countsValues[] = { id, PointerGetDatum(data[i].counts) };

            if(SPI_execute_plan(countsPlan, countsValues, NULL, false, 0) != SPI_OK_INSERT)
                elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

            for(int j = 0; j < COUNTS_SIZE; j++)
                counts[currentSeqid][j] = ((int16 *) ARR_DATA_PTR(data[i].counts))[j];
#endif

            Datum fingerprintValues[] = { id, Int16GetDatum(data[i].bitCount), PointerGetDatum(data[i].fp) };

            if(SPI_execute_plan(fingerprintPlan, fingerprintValues, NULL, false, 0) != SPI_OK_INSERT)
                elog(ERROR, "%s: SPI_execute_plan() failed", __func__);


            int length = ARR_DIMS(data[i].fp)[0];
            uint64_t *fpdata = (uint64_t *) ARR_DATA_PTR(data[i].fp);

            BitSet fp;
            bitset_init(&fp, fpdata, length);

            for(int j = bitset_next_set_bit(&fp, 1); j >= 0 && j <= FP_SIZE; j = bitset_next_set_bit(&fp, j + 1))
                bitset_set(bitmap + j - 1, currentSeqid);


            pfree(data[i].counts);
            pfree(data[i].fp);
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


    /*
     * create index
     */

    create_base_directory();

    if(unlikely(SPI_exec("select id from " INDEX_TABLE, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    int indexNumber = 0;

    if(SPI_processed != 0)
    {
        if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
            elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

        Datum number = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE) || isNullFlag)
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        indexNumber = DatumGetInt32(number) + 1;
    }


    char *indexFilePath = get_index_path(ORCHEM_INDEX_PREFIX, ORCHEM_INDEX_SUFFIX, indexNumber);

    if(unlikely(SPI_exec("delete from " INDEX_TABLE, 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(SPI_execute_with_args("insert into " INDEX_TABLE " (id) values ($1)", 1, (Oid[]) { INT4OID },
            (Datum[]) {Int32GetDatum(indexNumber)}, NULL, false, 0) != SPI_OK_INSERT)
        elog(ERROR, "%s: SPI_execute_with_args() failed", __func__);


    int fd = open(indexFilePath, O_EXCL | O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);

    if(fd == -1)
        elog(ERROR, "%s: open() failed", __func__);

    PG_TRY();
    {
        uint64_t offset = FP_SIZE + 1;
        int indexSize = maxIndexSize - notIndexed;

#if USE_COUNT_FINGERPRINT
        offset += (indexSize * COUNTS_SIZE * sizeof(uint16) + sizeof(uint64_t) - 1) / sizeof(uint64_t);
#endif

        for(int i = 0; i < FP_SIZE; i++)
        {
            if(write(fd, &offset, sizeof(uint64_t)) != sizeof(uint64_t))
                elog(ERROR, "%s: write() failed", __func__);

            offset += bitmap[i].wordsInUse + 1;
        }

        if(write(fd, &indexSize, sizeof(int64_t)) != sizeof(int64_t))
            elog(ERROR, "%s: write() failed", __func__);


#if USE_COUNT_FINGERPRINT
        for(int i = 0; i < COUNTS_SIZE; i++)
            if(write(fd, (int16 *) counts + i * indexSize, indexSize * sizeof(int16)) != (ssize_t) (indexSize * sizeof(int16)))
                elog(ERROR, "%s: write() failed", __func__);


        uint64_t zero = 0;
        ssize_t padding = (sizeof(uint64_t) - (indexSize * COUNTS_SIZE * sizeof(uint16)) % sizeof(uint64_t)) % sizeof(uint64_t);

        if(write(fd, &zero, padding) != padding)
            elog(ERROR, "%s: write() failed", __func__);
#endif


        for(int i = 0; i < FP_SIZE; i++)
        {
            uint64_t wordsInUse = bitmap[i].wordsInUse;

            if(write(fd, &wordsInUse, sizeof(uint64_t)) != sizeof(uint64_t))
                elog(ERROR, "%s: write() failed", __func__);

            if(write(fd, bitmap[i].words, wordsInUse * sizeof(uint64_t)) != (ssize_t) (wordsInUse * sizeof(uint64_t)))
                elog(ERROR, "%s: write() failed", __func__);
        }


        if(close(fd) != 0)
            elog(ERROR, "%s: close() failed", __func__);

#if USE_MOLECULE_INDEX
        sachem_generate_molecule_index(indexNumber, true);
#endif
    }
    PG_CATCH();
    {
        unlink(indexFilePath);

        PG_RE_THROW();
    }
    PG_END_TRY();


    SPI_finish();
    PG_RETURN_VOID();
}
