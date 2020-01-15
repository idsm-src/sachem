#include <postgres.h>
#include <catalog/pg_type.h>
#include <catalog/namespace.h>
#include <executor/spi.h>
#include <funcapi.h>
#include "enum.h"
#include "java.h"
#include "sachem.h"


#define BUFFERSIZE 1000000


typedef struct
{
    jobject handler;
    int32 limit;
    int32 size;
    int32 possition;
    int32 ids[BUFFERSIZE];
    float4 scores[BUFFERSIZE];
}
LuceneResult;


static bool initialized = false;
static TupleDesc tupdesc = NULL;
static SPIPlanPtr snapshotQueryPlan;

static EnumValue searchModeTable[2];
static EnumValue chargeModeTable[3];
static EnumValue isotopeModeTable[3];
static EnumValue stereoModeTable[2];
static EnumValue aromaticityModeTable[3];
static EnumValue tautomerModeTable[2];

static jmethodID setFolderMethod;
static jmethodID subsearchMethod;
static jmethodID simsearchMethod;
static jmethodID loadMethod;
static jfieldID idsField;
static jfieldID scoresField;

static jobject lucene = NULL;
static int indexId = -1;


static void lucene_set_folder(jobject lucene, const char *path)
{
    jstring folder = NULL;

    PG_TRY();
    {
        folder = (*env)->NewStringUTF(env, path);
        java_check_exception(__func__);

        (*env)->CallVoidMethod(env, lucene, setFolderMethod, folder);
        java_check_exception(__func__);

        JavaDeleteRef(folder);
    }
    PG_CATCH();
    {
        JavaDeleteRef(folder);

        PG_RE_THROW();
    }
    PG_END_TRY();
}


int lucene_search_update_snapshot(jobject lucene)
{
    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    if(unlikely(SPI_execute_plan(snapshotQueryPlan, NULL, NULL, true, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    bool isNullFlag;
    int32_t dbIndexNumber = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);

    SPI_finish();


    if(unlikely(dbIndexNumber != indexId))
    {
        char *path = get_index_path(INDEX_PREFIX_NAME, dbIndexNumber);
        lucene_set_folder(lucene, path);
        indexId = dbIndexNumber;
    }

    return indexId;
}


static void lucene_search_init(void)
{
    if(likely(initialized))
        return;

    java_init();

    Oid spaceid = LookupExplicitNamespace("sachem", false);


    /* search modes */
    {
        Oid typoid = LookupExplicitEnumType(spaceid, "search_mode");

        jclass enumClass = (*env)->FindClass(env, "cz/iocb/sachem/molecule/SearchMode");
        java_check_exception(__func__);

        const char* const values[] = { "SUBSTRUCTURE", "EXACT" };

        for(int i = 0; i < 2; i++)
        {
            searchModeTable[i].oid = LookupExplicitEnumValue(typoid, values[i]);
            searchModeTable[i].object = LookupJavaEnumValue(enumClass, values[i], "Lcz/iocb/sachem/molecule/SearchMode;");
        }
    }


    /* charge modes */
    {
        Oid typoid = LookupExplicitEnumType(spaceid, "charge_mode");

        jclass enumClass = (*env)->FindClass(env, "cz/iocb/sachem/molecule/ChargeMode");
        java_check_exception(__func__);

        const char* const values[] = { "IGNORE", "DEFAULT_AS_UNCHARGED", "DEFAULT_AS_ANY" };

        for(int i = 0; i < 3; i++)
        {
            chargeModeTable[i].oid = LookupExplicitEnumValue(typoid, values[i]);
            chargeModeTable[i].object = LookupJavaEnumValue(enumClass, values[i], "Lcz/iocb/sachem/molecule/ChargeMode;");
        }
    }


    /* isotope modes */
    {
        Oid typoid = LookupExplicitEnumType(spaceid, "isotope_mode");

        jclass enumClass = (*env)->FindClass(env, "cz/iocb/sachem/molecule/IsotopeMode");
        java_check_exception(__func__);

        const char* const values[] = { "IGNORE", "DEFAULT_AS_STANDARD", "DEFAULT_AS_ANY" };

        for(int i = 0; i < 3; i++)
        {
            isotopeModeTable[i].oid = LookupExplicitEnumValue(typoid, values[i]);
            isotopeModeTable[i].object = LookupJavaEnumValue(enumClass, values[i], "Lcz/iocb/sachem/molecule/IsotopeMode;");
        }
    }


    /* stereo modes */
    {
        Oid typoid = LookupExplicitEnumType(spaceid, "stereo_mode");

        jclass enumClass = (*env)->FindClass(env, "cz/iocb/sachem/molecule/StereoMode");
        java_check_exception(__func__);

        const char* const values[] = { "IGNORE", "STRICT" };

        for(int i = 0; i < 2; i++)
        {
            stereoModeTable[i].oid = LookupExplicitEnumValue(typoid, values[i]);
            stereoModeTable[i].object = LookupJavaEnumValue(enumClass, values[i], "Lcz/iocb/sachem/molecule/StereoMode;");
        }
    }


    /* aromaticity modes */
    {
        Oid typoid = LookupExplicitEnumType(spaceid, "aromaticity_mode");

        jclass enumClass = (*env)->FindClass(env, "cz/iocb/sachem/molecule/AromaticityMode");
        java_check_exception(__func__);

        const char* const values[] = { "PRESERVE", "DETECT", "AUTO" };

        for(int i = 0; i < 3; i++)
        {
            aromaticityModeTable[i].oid = LookupExplicitEnumValue(typoid, values[i]);
            aromaticityModeTable[i].object = LookupJavaEnumValue(enumClass, values[i], "Lcz/iocb/sachem/molecule/AromaticityMode;");
        }
    }


    /* tautomer modes */
    {
        Oid typoid = LookupExplicitEnumType(spaceid, "tautomer_mode");

        jclass enumClass = (*env)->FindClass(env, "cz/iocb/sachem/molecule/TautomerMode");
        java_check_exception(__func__);

        const char* const values[] = { "IGNORE", "INCHI" };

        for(int i = 0; i < 2; i++)
        {
            tautomerModeTable[i].oid = LookupExplicitEnumValue(typoid, values[i]);
            tautomerModeTable[i].object = LookupJavaEnumValue(enumClass, values[i], "Lcz/iocb/sachem/molecule/TautomerMode;");
        }
    }


    /* create tuple description */
    if(unlikely(tupdesc == NULL))
    {
        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
            elog(ERROR, "%s: SPI_connect() failed", __func__);

        TupleDesc desc = NULL;

        PG_MEMCONTEXT_BEGIN(TopMemoryContext);
        PG_TRY();
        {
            desc = CreateTemplateTupleDesc(2, false);
            TupleDescInitEntry(desc, (AttrNumber) 1, "compound", INT4OID, -1, 0);
            TupleDescInitEntry(desc, (AttrNumber) 2, "score", FLOAT4OID, -1, 0);
            desc = BlessTupleDesc(desc);
            tupdesc = desc;
        }
        PG_CATCH();
        {
            if(desc != NULL)
                FreeTupleDesc(desc);

            PG_RE_THROW();
        }
        PG_END_TRY();
        PG_MEMCONTEXT_END();

        SPI_finish();
    }


    /* prepare snapshot query plan */
    if(unlikely(snapshotQueryPlan == NULL))
    {
        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
            elog(ERROR, "%s: SPI_connect() failed", __func__);

        SPIPlanPtr plan = SPI_prepare("select id from " INDEX_TABLE, 0, NULL);

        if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        snapshotQueryPlan = plan;

        SPI_finish();
    }


    /* java init */
    jclass searcherClass = (*env)->FindClass(env, "cz/iocb/sachem/lucene/Searcher");
    java_check_exception(__func__);

    jmethodID constructor = (*env)->GetMethodID(env, searcherClass, "<init>", "()V");
    java_check_exception(__func__);

    setFolderMethod = (*env)->GetMethodID(env, searcherClass, "setFolder", "(Ljava/lang/String;)V");
    java_check_exception(__func__);

    subsearchMethod = (*env)->GetMethodID(env, searcherClass, "subsearch", "([BLcz/iocb/sachem/molecule/SearchMode;Lcz/iocb/sachem/molecule/ChargeMode;Lcz/iocb/sachem/molecule/IsotopeMode;Lcz/iocb/sachem/molecule/StereoMode;Lcz/iocb/sachem/molecule/AromaticityMode;Lcz/iocb/sachem/molecule/TautomerMode;I)Lcz/iocb/sachem/lucene/SearcherHandler;");
    java_check_exception(__func__);

    simsearchMethod = (*env)->GetMethodID(env, searcherClass, "simsearch", "([BFILcz/iocb/sachem/molecule/AromaticityMode;Lcz/iocb/sachem/molecule/TautomerMode;I)Lcz/iocb/sachem/lucene/SearcherHandler;");
    java_check_exception(__func__);


    jclass handlerClass = (*env)->FindClass(env, "cz/iocb/sachem/lucene/SearcherHandler");
    java_check_exception(__func__);

    loadMethod = (*env)->GetMethodID(env, handlerClass, "load", "(I)I");
    java_check_exception(__func__);

    idsField = (*env)->GetFieldID(env, handlerClass, "ids", "[I");
    java_check_exception(__func__);

    scoresField = (*env)->GetFieldID(env, handlerClass, "scores", "[F");
    java_check_exception(__func__);


    lucene = (*env)->NewObject(env, searcherClass, constructor);
    java_check_exception(__func__);


    initialized = true;
}


static HeapTuple lucene_result_get_item(LuceneResult *result)
{
    if(result->possition == BUFFERSIZE && result->limit != 0)
    {
        int limit = BUFFERSIZE;

        if(result->limit > 0)
        {
            if(result->limit < BUFFERSIZE)
                limit = result->limit;

            result->limit -= limit;
        }


        result->possition = 0;
        result->size = (*env)->CallIntMethod(env, result->handler, loadMethod, limit);
        java_check_exception(__func__);


        jarray array = NULL;

        PG_TRY();
        {
            array = (*env)->GetObjectField(env, result->handler, idsField);
            jint *ids = (*env)->GetPrimitiveArrayCritical(env, array, NULL);

            if(ids == NULL)
                elog(ERROR, "out of memmory");

            memcpy(result->ids, ids, result->size * sizeof(int32));

            (*env)->ReleasePrimitiveArrayCritical(env, array, ids, JNI_ABORT);
            JavaDeleteRef(array);


            array = (*env)->GetObjectField(env, result->handler, scoresField);
            jfloat *scores = (*env)->GetPrimitiveArrayCritical(env, array, NULL);

            if(scores == NULL)
                elog(ERROR, "out of memmory");

            memcpy(result->scores, scores, result->size * sizeof(float4));

            (*env)->ReleasePrimitiveArrayCritical(env, array, scores, JNI_ABORT);
            JavaDeleteRef(array);
        }
        PG_CATCH();
        {
            JavaDeleteRef(array);

            PG_RE_THROW();
        }
        PG_END_TRY();
    }


    if(result->possition == result->size)
        return NULL;


    bool isnull[2] = {0, 0};
    Datum values[2] = {Int32GetDatum(result->ids[result->possition]), Float4GetDatum(result->scores[result->possition])};
    HeapTuple tuple = heap_form_tuple(tupdesc, values, isnull);

    result->possition++;

    return tuple;
}


static void lucene_result_free(LuceneResult *result)
{
    if(result != NULL)
        JavaDeleteRef(result->handler);
}


static LuceneResult *lucene_subsearch(jobject lucene, VarChar *query, Oid search, Oid charge, Oid isotope, Oid stereo, Oid aromaticity, Oid tautomers, int32 limit)
{
    LuceneResult *result = NULL;
    jbyteArray queryArray = NULL;
    jobject handler = NULL;


    PG_TRY();
    {
        size_t length = VARSIZE(query) - VARHDRSZ;

        queryArray = (jbyteArray) (*env)->NewByteArray(env, length);
        java_check_exception(__func__);

        (*env)->SetByteArrayRegion(env, queryArray, 0, length, (jbyte *) VARDATA(query));
        java_check_exception(__func__);

        handler = (*env)->CallObjectMethod(env, lucene, subsearchMethod, queryArray,
                ConvertEnumValue(searchModeTable, search),
                ConvertEnumValue(chargeModeTable, charge),
                ConvertEnumValue(isotopeModeTable, isotope),
                ConvertEnumValue(stereoModeTable, stereo),
                ConvertEnumValue(aromaticityModeTable,  aromaticity),
                ConvertEnumValue(tautomerModeTable,  tautomers),
                BUFFERSIZE);
        java_check_exception(__func__);

        JavaDeleteRef(queryArray);

        result = palloc(sizeof(LuceneResult));
        result->handler = handler;
        result->limit = limit;
        result->size = BUFFERSIZE;
        result->possition = BUFFERSIZE;
    }
    PG_CATCH();
    {
        JavaDeleteRef(queryArray);
        JavaDeleteRef(handler);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return result;
}


static LuceneResult *lucene_simsearch(jobject lucene, VarChar *query, float4 threshold, int32 depth, Oid aromaticity, Oid tautomers, int32 limit)
{
    LuceneResult *result = NULL;
    jbyteArray queryArray = NULL;
    jobject handler = NULL;


    PG_TRY();
    {
        size_t length = VARSIZE(query) - VARHDRSZ;

        queryArray = (jbyteArray) (*env)->NewByteArray(env, length);
        java_check_exception(__func__);

        (*env)->SetByteArrayRegion(env, queryArray, 0, length, (jbyte *) VARDATA(query));
        java_check_exception(__func__);


        handler = (*env)->CallObjectMethod(env, lucene, simsearchMethod, queryArray, threshold, depth,
                ConvertEnumValue(aromaticityModeTable, aromaticity),
                ConvertEnumValue(tautomerModeTable, tautomers),
                BUFFERSIZE);
        java_check_exception(__func__);

        JavaDeleteRef(queryArray);

        result = palloc(sizeof(LuceneResult));
        result->handler = handler;
        result->limit = limit;
        result->size = BUFFERSIZE;
        result->possition = BUFFERSIZE;
    }
    PG_CATCH();
    {
        JavaDeleteRef(queryArray);
        JavaDeleteRef(handler);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return result;
}


PG_FUNCTION_INFO_V1(sachem_substructure_search);
Datum sachem_substructure_search(PG_FUNCTION_ARGS)
{
    if(SRF_IS_FIRSTCALL())
    {
        lucene_search_init();
        lucene_search_update_snapshot(lucene);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        Oid search = PG_GETARG_OID(1);
        Oid charge = PG_GETARG_OID(2);
        Oid isotope = PG_GETARG_OID(3);
        Oid stereo = PG_GETARG_OID(4);
        Oid aromaticity = PG_GETARG_OID(5);
        Oid tautomers = PG_GETARG_OID(6);
        int32 limit = PG_GETARG_INT32(7);

        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);
        funcctx->user_fctx = lucene_subsearch(lucene, query, search, charge, isotope, stereo, aromaticity, tautomers, limit);
        PG_MEMCONTEXT_END();
    }


    FuncCallContext *funcctx = SRF_PERCALL_SETUP();
    LuceneResult *result = funcctx->user_fctx;
    HeapTuple item;

    PG_TRY();
    {
        item = lucene_result_get_item(result);
    }
    PG_CATCH();
    {
        lucene_result_free(result);
        PG_RE_THROW();
    }
    PG_END_TRY();

    if(HeapTupleIsValid(item))
        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(item));

    lucene_result_free(result);
    SRF_RETURN_DONE(funcctx);
}


PG_FUNCTION_INFO_V1(sachem_similarity_search);
Datum sachem_similarity_search(PG_FUNCTION_ARGS)
{
    if(SRF_IS_FIRSTCALL())
    {
        lucene_search_init();
        lucene_search_update_snapshot(lucene);

        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();

        VarChar *query = PG_GETARG_VARCHAR_P(0);
        float4 threshold = PG_GETARG_FLOAT4(1);
        int32 depth = PG_GETARG_INT32(2);
        Oid aromaticity = PG_GETARG_OID(3);
        Oid tautomers = PG_GETARG_OID(4);
        int32 limit = PG_GETARG_INT32(5);

        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);
        funcctx->user_fctx = lucene_simsearch(lucene, query, threshold, depth, aromaticity, tautomers, limit);
        PG_MEMCONTEXT_END();
    }


    FuncCallContext *funcctx = SRF_PERCALL_SETUP();
    LuceneResult *result = funcctx->user_fctx;
    HeapTuple item;

    PG_TRY();
    {
        item = lucene_result_get_item(result);
    }
    PG_CATCH();
    {
        lucene_result_free(result);
        PG_RE_THROW();
    }
    PG_END_TRY();

    if(HeapTupleIsValid(item))
        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(item));

    lucene_result_free(result);
    SRF_RETURN_DONE(funcctx);
}
