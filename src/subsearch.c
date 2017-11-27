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
#define COUNTS_SIZE               10
#define MOLECULES_TABLE           "orchem_molecules"
#define MOLECULE_COUNTS_TABLE     "orchem_molecule_counts"
#define FINGERPRINT_INDEX_TABLE   "orchem_substructure_fingerprint_index"
#define COUNT_INDEX_FILE          "orchem_molecule_counts.idx"
#define FINGERPRINT_INDEX_FILE    "orchem_substructure_fingerprint_index.idx"


typedef struct
{
    int32_t topN;
    bool strictStereo;
    bool exact;
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


static bool initialised = false;
static int moleculeCount;
static BitSet bitmap[FP_SIZE];
static SPIPlanPtr mainQueryPlan;

#if USE_COUNT_FINGERPRINT
static int16 *counts[COUNTS_SIZE];
#endif


void subsearch_module_init(void)
{
    /* load bit fingerprint */
    char *fingerprintIndexFilePath = getFilePath(FINGERPRINT_INDEX_FILE);

    if(access(fingerprintIndexFilePath, R_OK) != -1)
    {
        int fd = open(fingerprintIndexFilePath, O_RDONLY, 0);

        struct stat st;
        fstat(fd, &st);
        size_t length = st.st_size;

        uint64_t *addr = mmap(NULL, length, PROT_READ, MAP_SHARED, fd, 0);
        moleculeCount = *(addr++);

        for(int i = 0; i < FP_SIZE; i++)
        {
            int length = *addr;

            bitset_init(bitmap + i, addr + 1, length);
            addr += length + 1;
        }
    }


#if USE_COUNT_FINGERPRINT
    /* load count fingerprint */
    char *countIndexFilePath = getFilePath(COUNT_INDEX_FILE);

    if(access(countIndexFilePath, R_OK) != -1)
    {
        int fd = open(countIndexFilePath, O_RDONLY, 0);

        struct stat st;
        fstat(fd, &st);
        size_t length = st.st_size;

        void *addr = mmap(NULL, length, PROT_READ, MAP_SHARED, fd, 0);

        if(*((uint64_t *) addr) != moleculeCount)
            elog(ERROR, "subsearch module: wrong index size");

        addr += sizeof(uint64_t);

        for(int j = 0; j < COUNTS_SIZE; j++)
            counts[j] =  ((int16 *) addr) + j * moleculeCount;
    }
#endif


    PG_TRY();
    {
        /* prepare query plan */
        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
            elog(ERROR, "subsearch module: SPI_connect() failed");

        mainQueryPlan = SPI_prepare("select id, seqid, atoms, bonds from " MOLECULES_TABLE " where seqid = any($1)", 1, (Oid[]) { INT4ARRAYOID });

        if(unlikely(mainQueryPlan == NULL))
            elog(ERROR, "subsearch module: SPI_prepare_cursor() failed");

        if(unlikely(SPI_keepplan(mainQueryPlan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "subsearch module: SPI_keepplan() failed");


        SPI_finish();
        initialised = true;
    }
    PG_CATCH();
    {
        elog(NOTICE, "subsearch module: initialization failed");
    }
    PG_END_TRY();
}


void subsearch_module_finish(void)
{
    initialised = false;
}


PG_FUNCTION_INFO_V1(orchem_substructure_search);
Datum orchem_substructure_search(PG_FUNCTION_ARGS)
{
    if(SRF_IS_FIRSTCALL())
    {
#if SHOW_STATS
        struct timeval begin;
        gettimeofday(&begin, NULL);
#endif

        if(unlikely(!initialised))
            elog(ERROR, "subsearch module is not properly initialized");

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        text *type = PG_GETARG_TEXT_P(1);
        int32_t topN = PG_GETARG_INT32(2);
        bool strictStereo = PG_GETARG_BOOL(3);
        bool exact = PG_GETARG_BOOL(4);
        bool tautomers = PG_GETARG_BOOL(5);
        int32_t vf2_timeout = PG_GETARG_INT32(6);
        char *typeStr = text_to_cstring(type);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();
        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);

        SubstructureSearchData *info = (SubstructureSearchData *) palloc(sizeof(SubstructureSearchData));
        funcctx->user_fctx = info;

        info->topN = topN;
        info->strictStereo = strictStereo;
        info->exact = exact;
        info->vf2_timeout = vf2_timeout;

        info->queryDataCount = java_parse_orchem_substructure_query(&info->queryData, VARDATA(query), VARSIZE(query) - VARHDRSZ, typeStr, tautomers);

        PG_FREE_IF_COPY(query, 0);
        PG_FREE_IF_COPY(type, 1);
        pfree(typeStr);

        info->queryDataPosition = -1;
        info->candidatePosition = -1;
        info->table = NULL;
        info->tableRowCount = -1;
        info->tableRowPosition = -1;
        info->foundResults = 0;

        bitset_init_alloc(&info->candidates, moleculeCount);
        bitset_init_setted(&info->resultMask, moleculeCount);

        info->isomorphismContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx, "subsearch isomorphism context",
                ALLOCSET_DEFAULT_MINSIZE, ALLOCSET_DEFAULT_INITSIZE, ALLOCSET_DEFAULT_MAXSIZE);
        info->targetContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx, "subsearch target context",
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
                    molecule_init(&info->queryMolecule, data->atomLength, data->atoms, data->bondLength, data->bonds, data->restH, !info->exact);
                    vf2state_init(&info->vf2state, &info->queryMolecule, info->strictStereo, info->exact);
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
                        if(counts[j][info->candidatePosition] < data->counts[j])
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
                     elog(ERROR, "subsearch module: SPI_connect() failed");

                connected = true;


                if(unlikely(SPI_execute_plan(mainQueryPlan, values, NULL, true, 0) != SPI_OK_SELECT))
                    elog(ERROR, "subsearch module: SPI_execute_plan() failed");

                if(unlikely(SPI_processed != count || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 4))
                    elog(ERROR, "subsearch module: SPI_execute_plan() failed");

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
                elog(ERROR, "subsearch module: SPI_getbinval() failed");


            Datum seqid = SPI_getbinval(tuple, tupdesc, 2, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch module: SPI_getbinval() failed");


            Datum atoms = SPI_getbinval(tuple, tupdesc, 3, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch module: SPI_getbinval() failed");


            Datum bonds = SPI_getbinval(tuple, tupdesc, 4, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch module: SPI_getbinval() failed");

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


#if USE_COUNT_FINGERPRINT
void orchem_substructure_write_count_index()
{
    char *indexFilePath = getFilePath(COUNT_INDEX_FILE);
    int fd = open(indexFilePath, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);

    if(fd == -1)
        elog(ERROR, "orchem_substructure_write_count_index: open() failed");


    char isNullFlag;
    MemoryContext context = CurrentMemoryContext;

    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "orchem_substructure_write_count_index: SPI_connect() failed");


    if(unlikely(SPI_execute("select max(seqid)+1 from " MOLECULE_COUNTS_TABLE, true, FETCH_ALL) != SPI_OK_SELECT))
        elog(ERROR, "orchem_substructure_write_count_index: SPI_execute() failed");

    if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
        elog(ERROR, "orchem_substructure_write_count_index: SPI_execute() failed");

    int64_t moleculeCount = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "orchem_substructure_write_count_index: SPI_getbinval() failed");

    if(write(fd, &moleculeCount, sizeof(int64_t)) != sizeof(int64_t))
        elog(ERROR, "orchem_substructure_write_count_index: write() failed");


    int16 *counts[COUNTS_SIZE];
    PG_MEMCONTEXT_BEGIN(context);

    for(int j = 0; j < COUNTS_SIZE; j++)
        counts[j] = (int16 *) palloc0(moleculeCount * sizeof(int16));

    PG_MEMCONTEXT_END();


    Portal cursor = SPI_cursor_open_with_args(NULL, "select seqid, molTripleBondCount, molSCount, molOCount, molNCount, molFCount, molClCount, molBrCount, "
            "molICount, molCCount, molPCount from " MOLECULE_COUNTS_TABLE " order by seqid",
            0, NULL, NULL, NULL, true, CURSOR_OPT_BINARY | CURSOR_OPT_NO_SCROLL);

    if(unlikely(cursor == NULL))
            elog(ERROR, "orchem_substructure_write_count_index: SPI_cursor_open_with_args() failed");

    while(true)
    {
        SPI_cursor_fetch(cursor, true, 100000);

        if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != COUNTS_SIZE + 1))
            elog(ERROR, "orchem_substructure_write_count_index: SPI_cursor_fetch() failed");

        if(SPI_processed == 0)
            break;

        for(int i = 0; i < SPI_processed; i++)
        {
            HeapTuple tuple = SPI_tuptable->vals[i];

            int32 seqid = DatumGetInt32(SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag));

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "orchem_substructure_write_count_index: SPI_getbinval() failed");

            for(int j = 0; j < COUNTS_SIZE; j++)
            {
                int16 value = DatumGetInt16(SPI_getbinval(tuple, SPI_tuptable->tupdesc, j + 2, &isNullFlag));

                if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                    elog(ERROR, "orchem_substructure_write_count_index: SPI_getbinval() failed");

                counts[j][seqid] = value;
            }
        }

        SPI_freetuptable(SPI_tuptable);
    }


    for(int j = 0; j < COUNTS_SIZE; j++)
    {
        ssize_t writen = write(fd, counts[j], moleculeCount * sizeof(int16));

        if(writen != moleculeCount * sizeof(int16))
            elog(ERROR, "orchem_substructure_write_count_index: write() failed");
    }


    SPI_cursor_close(cursor);
    SPI_finish();

    if(close(fd) != 0)
        elog(ERROR, "orchem_substructure_write_count_index: close() failed");
}
#endif


void orchem_substructure_write_fp_index()
{
    char *indexFilePath = getFilePath(FINGERPRINT_INDEX_FILE);
    int fd = open(indexFilePath, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);

    if(fd == -1)
        elog(ERROR, "orchem_substructure_write_fp_index: open() failed");


    char isNullFlag;

    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "orchem_substructure_write_fp_index: SPI_connect() failed");


    if(unlikely(SPI_execute("select max(seqid)+1 from " MOLECULE_COUNTS_TABLE, true, FETCH_ALL) != SPI_OK_SELECT))
        elog(ERROR, "orchem_substructure_write_fp_index: SPI_execute() failed");

    if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
        elog(ERROR, "orchem_substructure_write_fp_index: SPI_execute() failed");

    int64_t moleculeCount = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "orchem_substructure_write_fp_index: SPI_getbinval() failed");

    if(write(fd, &moleculeCount, sizeof(int64_t)) != sizeof(int64_t))
        elog(ERROR, "orchem_substructure_write_fp_index: write() failed");


    if(unlikely(SPI_execute("select bitmap from " FINGERPRINT_INDEX_TABLE " order by idx", true, FETCH_ALL) != SPI_OK_SELECT))
        elog(ERROR, "orchem_substructure_write_fp_index: SPI_execute() failed");

    if(unlikely(SPI_processed != FP_SIZE || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
        elog(ERROR, "orchem_substructure_write_fp_index: SPI_execute() failed");

    for(int i = 0; i < SPI_processed; i++)
    {
        HeapTuple tuple = SPI_tuptable->vals[i];

        bytea *fp = DatumGetByteaP(SPI_getbinval(tuple, SPI_tuptable->tupdesc, 1, &isNullFlag));

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "orchem_substructure_write_fp_index: SPI_getbinval() failed");

        BitSet bitmap;
        bitset_init_from_array(&bitmap, VARDATA(fp), VARSIZE(fp) - VARHDRSZ);

        uint64_t wordsInUse = bitmap.wordsInUse;

        if(write(fd, &wordsInUse, sizeof(uint64_t)) != sizeof(uint64_t))
            elog(ERROR, "orchem_substructure_write_count_index: write() failed");

        if(write(fd, bitmap.words, wordsInUse * sizeof(uint64_t)) != wordsInUse * sizeof(uint64_t))
            elog(ERROR, "orchem_substructure_write_count_index: write() failed");
    }


    SPI_freetuptable(SPI_tuptable);
    SPI_finish();

    if(close(fd) != 0)
        elog(ERROR, "orchem_substructure_write_fp_index: close() failed");
}


PG_FUNCTION_INFO_V1(orchem_substructure_write_indexes);
Datum orchem_substructure_write_indexes(PG_FUNCTION_ARGS)
{
    createBasePath();

    orchem_substructure_write_fp_index();

#if USE_COUNT_FINGERPRINT
    orchem_substructure_write_count_index();
#endif

    PG_RETURN_BOOL(true);
}
