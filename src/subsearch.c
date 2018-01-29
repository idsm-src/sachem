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
#include <unistd.h>
#include "bitset.h"
#include "isomorphism.h"
#include "java.h"
#include "molecule.h"
#include "sachem.h"
#include "subsearch.h"


#define SHOW_STATS                0
#define USE_COUNT_FINGERPRINT     1
#define FETCH_SIZE                1000
#define FP_SIZE                   788
#define COUNTS_SIZE               13
#define SYNC_FETCH_SIZE           100000
#define COMPOUNDS_TABLE           "compounds"
#define AUDIT_TABLE               "orchem_compound_audit"
#define INDEX_TABLE               "orchem_index"
#define MOLECULES_TABLE           "orchem_molecules"
#define MOLECULE_ERRORS_TABLE     "orchem_molecule_errors"
#define FINGERPRINT_TABLE         "orchem_fingerprint"
#define MOLECULE_COUNTS_TABLE     "orchem_molecule_counts"


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

    OrchemSubstructureQueryData *queryData;
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


static int indexId = -1;
static uint64_t *base = MAP_FAILED;
static size_t length;
static int moleculeCount;
static BitSet bitmap[FP_SIZE];
static SPIPlanPtr indexQueryPlan = NULL;
static SPIPlanPtr mainQueryPlan = NULL;

#if USE_COUNT_FINGERPRINT
static int16 (*counts)[COUNTS_SIZE];
#endif


static void orchem_subsearch_init(void)
{
    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);


    /* prepare index query plan */
    if(unlikely(indexQueryPlan == NULL))
    {
        SPIPlanPtr plan = SPI_prepare("select id, path from " INDEX_TABLE, 0, NULL);

        if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        indexQueryPlan = plan;
    }


    /* main query plan */
    if(unlikely(mainQueryPlan == NULL))
    {
        SPIPlanPtr plan = SPI_prepare("select id, seqid, molecule from " MOLECULES_TABLE " where seqid = any($1)", 1, (Oid[]) { INT4ARRAYOID });

        if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        mainQueryPlan = plan;
    }


    /* get index information */
    if(unlikely(SPI_execute_plan(indexQueryPlan, NULL, NULL, true, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);


    char isNullFlag;

    Datum id = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);

    if(unlikely(DatumGetInt32(id) != indexId))
    {
        if(likely(base != MAP_FAILED))
        {
            if(unlikely(munmap(base, length) < 0))
                elog(ERROR, "%s: munmap() failed", __func__);

            base = MAP_FAILED;
        }


        Datum path = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isNullFlag);

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "%s: SPI_getbinval() failed", __func__);

        int fd = -1;
        void *address = MAP_FAILED;
        struct stat st;

        PG_TRY();
        {
            if(unlikely((fd = open(TextDatumGetCString(path), O_RDONLY, 0)) < 0))
                elog(ERROR, "%s: open() failed", __func__);

            if(fstat(fd, &st) < 0)
                elog(ERROR, "%s: fstat() failed", __func__);

            if(unlikely((address = mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED))
                elog(ERROR, "%s: mmap() failed", __func__);

            if(unlikely(close(fd) < 0))
                elog(ERROR, "%s: close() failed", __func__);
        }
        PG_CATCH();
        {
            if(address != MAP_FAILED)
                munmap(address, st.st_size);

            if(fd != -1)
                close(fd);

            PG_RE_THROW();
        }
        PG_END_TRY();

        base = address;
        length = st.st_size;
        moleculeCount = *(base + FP_SIZE);

        for(int i = 0; i < FP_SIZE; i++)
        {
            uint64_t *address = base + base[i];
            bitset_init(bitmap + i, address + 1, *address);
        }

#if USE_COUNT_FINGERPRINT
        counts = (int16 (*)[COUNTS_SIZE]) (base + FP_SIZE + 1);
#endif

        indexId = DatumGetInt32(id);
    }

    SPI_finish();
}


PG_FUNCTION_INFO_V1(orchem_substructure_search);
Datum orchem_substructure_search(PG_FUNCTION_ARGS)
{
    if(SRF_IS_FIRSTCALL())
    {
        orchem_subsearch_init();

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
        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);

        SubstructureSearchData *info = (SubstructureSearchData *) palloc(sizeof(SubstructureSearchData));
        funcctx->user_fctx = info;

        info->topN = topN;
        info->graphMode = graphMode;
        info->chargeMode = chargeMode;
        info->isotopeMode = isotopeMode;
        info->stereoMode = stereoMode;
        info->vf2_timeout = vf2_timeout;

        info->queryDataCount = java_parse_orchem_substructure_query(&info->queryData, VARDATA(query), VARSIZE(query) - VARHDRSZ,
                type, graphMode == GRAPH_EXACT, tautomerMode == TAUTOMER_INCHI);

        PG_FREE_IF_COPY(query, 0);

        info->queryDataPosition = -1;
        info->candidatePosition = -1;
        info->table = NULL;
        info->tableRowCount = -1;
        info->tableRowPosition = -1;
        info->foundResults = 0;

        bitset_init_alloc(&info->candidates, moleculeCount);
        bitset_init_setted(&info->resultMask, moleculeCount);

        info->isomorphismContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx,
                "subsearch isomorphism context", ALLOCSET_DEFAULT_SIZES);
        info->targetContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx,
                "subsearch target context", ALLOCSET_DEFAULT_SIZES);

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
                {
                    MemoryContextDelete(info->table->tuptabcxt);
                    info->table = NULL;
                }

                if(info->candidatePosition < 0)
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

                    OrchemSubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                    MemoryContextReset(info->isomorphismContext);

                    PG_MEMCONTEXT_BEGIN(info->isomorphismContext);
                    info->extended = molecule_is_extended_search_needed(data->molecule, info->chargeMode, info->isotopeMode);
                    molecule_init(&info->queryMolecule, data->molecule, data->restH, info->extended,
                            info->chargeMode != CHARGE_IGNORE, info->isotopeMode != ISOTOPE_IGNORE, info->stereoMode != STEREO_IGNORE);
                    vf2state_init(&info->vf2state, &info->queryMolecule, info->graphMode, info->chargeMode,
                            info->isotopeMode, info->stereoMode);
                    PG_MEMCONTEXT_END();

                    info->candidatePosition = bitset_next_set_bit(&info->candidates, 0);
                }


                int32 *arrayData = (int32 *) ARR_DATA_PTR(info->arrayBuffer);
                OrchemSubstructureQueryData *data = &(info->queryData[info->queryDataPosition]);

                int count = 0;

                while(count < FETCH_SIZE && info->candidatePosition >= 0)
                {
#if USE_COUNT_FINGERPRINT
                    bool isValid = true;

                    for(int j = 0; j < COUNTS_SIZE; j++)
                    {
                        if(counts[info->candidatePosition][j] < data->counts[j])
                        {
                            isValid = false;
                            break;
                        }
                    }

                    if(isValid)
#endif
                        arrayData[count++] = info->candidatePosition;

                    info->candidatePosition = bitset_next_set_bit(&info->candidates, info->candidatePosition + 1);
                }

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

                if(unlikely(/*SPI_processed != count ||*/ SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 3))
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


            Datum seqid = SPI_getbinval(tuple, tupdesc, 2, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "%s: SPI_getbinval() failed", __func__);


            Datum molecule = SPI_getbinval(tuple, tupdesc, 3, &isNullFlag);

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
                bitset_unset(&info->resultMask, DatumGetInt32(seqid));
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


PG_FUNCTION_INFO_V1(orchem_sync_data);
Datum orchem_sync_data(PG_FUNCTION_ARGS)
{
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


        for(int i = 0; i < SPI_processed; i++)
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


        for(int i = 0; i < SPI_processed; i++)
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


        for(int i = 0; i < SPI_processed; i++)
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


        java_parse_orchem_data(processed, molfiles, data);


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
     *  create index
     */

    createBasePath();

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

    Name database = DatumGetName(DirectFunctionCall1(current_database, 0));
    size_t basePathLength = strlen(DataDir);
    size_t databaseLength = strlen(database->data);

    char *indexFilePath = (char *) palloc(basePathLength +  databaseLength + 64);
    sprintf(indexFilePath, "%s/%s/orchem_substructure_index-%i.idx", DataDir, database->data, indexNumber);


    if(unlikely(SPI_exec("delete from " INDEX_TABLE, 0) != SPI_OK_DELETE))
        elog(ERROR, "%s: SPI_exec() failed", __func__);

    if(SPI_execute_with_args("insert into " INDEX_TABLE " (id, path) values ($1,$2)", 2, (Oid[]) { INT4OID, TEXTOID },
            (Datum[]) {Int32GetDatum(indexNumber), CStringGetTextDatum(indexFilePath)}, NULL, false, 0) != SPI_OK_INSERT)
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
            if(write(fd, (int16 *) counts + i * indexSize, indexSize * sizeof(int16)) != indexSize * sizeof(int16))
                elog(ERROR, "%s: write() failed", __func__);


        uint64_t zero = 0;
        size_t padding = (sizeof(uint64_t) - (indexSize * COUNTS_SIZE * sizeof(uint16)) % sizeof(uint64_t)) % sizeof(uint64_t);

        if(write(fd, &zero, padding) != padding)
            elog(ERROR, "%s: write() failed", __func__);
#endif


        for(int i = 0; i < FP_SIZE; i++)
        {
            uint64_t wordsInUse = bitmap[i].wordsInUse;

            if(write(fd, &wordsInUse, sizeof(uint64_t)) != sizeof(uint64_t))
                elog(ERROR, "%s: write() failed", __func__);

            if(write(fd, bitmap[i].words, wordsInUse * sizeof(uint64_t)) != wordsInUse * sizeof(uint64_t))
                elog(ERROR, "%s: write() failed", __func__);
        }


        if(close(fd) != 0)
            elog(ERROR, "%s: close() failed", __func__);
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
