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
#include "subsearch-gin.h"


#define SHOW_STATS              0
#define FETCH_SIZE              1000
#define MOLECULES_TABLE         "orchem_molecules"
#define MOLECULE_COUNTS_TABLE   "orchem_molecule_counts"
#define FINGERPRINT_TABLE       "orchem_substructure_fingerprint"


typedef struct
{
    int32_t topN;
    bool strictStereo;
    bool exact;

    BitSet resultMask;
    int32_t foundResults;

    QueryData *queryData;
    int queryDataCount;
    int queryDataPosition;

    char *cursorName;

    SPITupleTable *table;
    int tableRowCount;
    int tableRowPosition;

    Molecule queryMolecule;
    VF2State vf2state;

    MemoryContext isomorphismContext;
    MemoryContext targetContext;

#if SHOW_STATS
    int candidateCount;
    struct timeval begin;
#endif

} SubstructureSearchData;


static bool initialised = false;
static SPIPlanPtr mainQueryPlan;


void subsearch_gin_module_init(void)
{
    PG_TRY();
    {
        /* prepare query plan */
        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
            elog(ERROR, "subsearch-gin module: SPI_connect() failed");

        mainQueryPlan = SPI_prepare("select m.id, m.seqid, m.atoms, m.bonds from " MOLECULES_TABLE " m, "
                FINGERPRINT_TABLE " f where m.seqid = f.seqid and f.fp @> $1", 1, (Oid[]) { INT2ARRAYOID });

        if(unlikely(mainQueryPlan == NULL))
            elog(ERROR, "subsearch-gin module: SPI_prepare_cursor() failed");

        if(unlikely(SPI_keepplan(mainQueryPlan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "subsearch-gin module: SPI_keepplan() failed");


        SPI_finish();
        initialised = true;
    }
    PG_CATCH();
    {
        elog(NOTICE, "subsearch-gin module: initialization failed");
    }
    PG_END_TRY();
}


void subsearch_gin_module_finish(void)
{
    initialised = false;
}


PG_FUNCTION_INFO_V1(orchem_substructure_gin_search);
Datum orchem_substructure_gin_search(PG_FUNCTION_ARGS)
{
    bool connected = false;


    if(SRF_IS_FIRSTCALL())
    {
#if SHOW_STATS
        struct timeval begin;
        gettimeofday(&begin, NULL);
#endif

        if(unlikely(!initialised))
            elog(ERROR, "subsearch-gin module is not properly initialized");

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        text *type = PG_GETARG_TEXT_P(1);
        int32_t topN = PG_GETARG_INT32(2);
        bool strictStereo = PG_GETARG_BOOL(3);
        bool exact = PG_GETARG_BOOL(4);
        bool tautomers = PG_GETARG_BOOL(5);
        char *typeStr = text_to_cstring(type);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();


        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
             elog(ERROR, "subsearch-gin module: SPI_connect() failed");

        connected = true;

        if(unlikely(SPI_execute("select max(seqid) + 1 from " FINGERPRINT_TABLE, true, FETCH_ALL) != SPI_OK_SELECT))
            elog(ERROR, "subsearch-gin module: SPI_execute() failed");

        if(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1)
            elog(ERROR, "subsearch-gin module: SPI_execute() failed");

        char isNullFlag;
        int64_t moleculeCount = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag));

        if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
            elog(ERROR, "subsearch-gin module: SPI_getbinval() failed");

        SPI_freetuptable(SPI_tuptable);


        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);

        SubstructureSearchData *info = palloc(sizeof(SubstructureSearchData));
        funcctx->user_fctx = info;

        info->topN = topN;
        info->strictStereo = strictStereo;
        info->exact = exact;

        info->queryDataCount = java_parse_query(&info->queryData, VARDATA(query), VARSIZE(query) - VARHDRSZ, typeStr, tautomers);

        PG_FREE_IF_COPY(query, 0);
        PG_FREE_IF_COPY(type, 1);
        pfree(typeStr);

        info->queryDataPosition = -1;
        info->cursorName = NULL;
        info->table = NULL;
        info->tableRowCount = -1;
        info->tableRowPosition = -1;
        info->foundResults = 0;

        bitset_init_empty(&info->resultMask, moleculeCount);

        info->isomorphismContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx, "subsearch-gin isomorphism context",
                ALLOCSET_DEFAULT_MINSIZE, ALLOCSET_DEFAULT_INITSIZE, ALLOCSET_DEFAULT_MAXSIZE);
        info->targetContext = AllocSetContextCreate(funcctx->multi_call_memory_ctx, "subsearch-gin target context",
                ALLOCSET_DEFAULT_MINSIZE, ALLOCSET_DEFAULT_INITSIZE, ALLOCSET_DEFAULT_MAXSIZE);

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
    Portal cursor = NULL;

    if(likely(info->topN <= 0 || info->topN != info->foundResults))
    {
        while(true)
        {
            if(info->table == NULL)
            {
                info->queryDataPosition++;

                if(unlikely(info->queryDataPosition == info->queryDataCount))
                    break;


                QueryData *data = &info->queryData[info->queryDataPosition];

                MemoryContextReset(info->isomorphismContext);

                PG_MEMCONTEXT_BEGIN(info->isomorphismContext);
                molecule_init(&info->queryMolecule, data->atomLength, data->atoms, data->bondLength, data->bonds, data->restH, !info->exact);
                vf2state_init(&info->vf2state, &info->queryMolecule, info->strictStereo, info->exact);
                PG_MEMCONTEXT_END();


                ArrayType *array = (ArrayType *) palloc(data->fpLength * sizeof(int16) + ARR_OVERHEAD_NONULLS(1));
                array->ndim = 1;
                array->dataoffset = 0;
                array->elemtype = INT2OID;
                memcpy(ARR_DATA_PTR(array), data->fp, data->fpLength * sizeof(int16));
                *(ARR_DIMS(array)) = data->fpLength;
                SET_VARSIZE(array, data->fpLength * sizeof(int16) + ARR_OVERHEAD_NONULLS(1));

                Datum values[] = { PointerGetDatum(array)};


                if(unlikely(!connected && SPI_connect() != SPI_OK_CONNECT))
                     elog(ERROR, "subsearch-gin module: SPI_connect() failed");

                connected = true;

                cursor = SPI_cursor_open(NULL, mainQueryPlan, values, NULL, true);

                if(unlikely(cursor == NULL))
                    elog(ERROR, "subsearch-gin module: SPI_cursor_open() failed");


                PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);
                size_t nameLength = strlen(cursor->name) + 1;
                info->cursorName = palloc(nameLength);
                memcpy(info->cursorName, cursor->name, nameLength);
                PG_MEMCONTEXT_END();

                info->tableRowPosition = 0;
                info->tableRowCount = 0;
            }

            if(unlikely(info->tableRowPosition == info->tableRowCount))
            {
                if(info->table != NULL)
                    MemoryContextDelete(info->table->tuptabcxt);


                if(unlikely(!connected && SPI_connect() != SPI_OK_CONNECT))
                     elog(ERROR, "subsearch-gin module: SPI_connect() failed");

                connected = true;

                if(cursor == NULL)
                {
                    cursor = SPI_cursor_find(info->cursorName);

                    if(unlikely(cursor == NULL))
                        elog(ERROR, "subsearch-gin module: SPI_cursor_find() failed");
                }

                SPI_cursor_fetch(cursor, true, FETCH_SIZE);

                if(unlikely(SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 4))
                    elog(ERROR, "subsearch-gin module: SPI_cursor_fetch() failed");

                info->table = SPI_tuptable;
                info->tableRowPosition = 0;
                info->tableRowCount = SPI_processed;

                MemoryContextSetParent(SPI_tuptable->tuptabcxt, funcctx->multi_call_memory_ctx);


                if(SPI_processed == 0)
                {
                    SPI_cursor_close(cursor);
                    MemoryContextDelete(info->table->tuptabcxt);

                    cursor = NULL;
                    info->table = NULL;
                    continue;
                }
            }


            TupleDesc tupdesc = info->table->tupdesc;
            HeapTuple tuple = info->table->vals[info->tableRowPosition++];
            char isNullFlag;

            Datum id = SPI_getbinval(tuple, tupdesc, 1, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch-gin module: SPI_getbinval() failed");


            int32 seqid = DatumGetInt32(SPI_getbinval(tuple, tupdesc, 2, &isNullFlag));

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch-gin module: SPI_getbinval() failed");


            if(bitset_get(&info->resultMask, seqid))
                continue;

            Datum atoms = SPI_getbinval(tuple, tupdesc, 3, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch-gin module: SPI_getbinval() failed");


            Datum bonds = SPI_getbinval(tuple, tupdesc, 4, &isNullFlag);

            if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
                elog(ERROR, "subsearch-gin module: SPI_getbinval() failed");

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
