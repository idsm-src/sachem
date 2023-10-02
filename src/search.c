#include <postgres.h>
#include <catalog/pg_type.h>
#include <catalog/namespace.h>
#include <executor/spi.h>
#include <utils/memutils.h>
#include <funcapi.h>
#include <math.h>
#include "enum.h"
#include "java.h"
#include "sachem.h"


typedef struct
{
    VarChar *index;
    jstring name;

    int32 length;
    int32 possition;

    jarray idsArray;
    jint *ids;

    jarray scoresArray;
    jfloat *scores;
}
LuceneResult;


static bool initialized = false;
static TupleDesc tupdesc = NULL;
static SPIPlanPtr configQueryPlan = NULL;

static EnumValue searchModeTable[2];
static EnumValue chargeModeTable[3];
static EnumValue isotopeModeTable[3];
static EnumValue radicalModeTable[3];
static EnumValue stereoModeTable[2];
static EnumValue aromaticityModeTable[3];
static EnumValue tautomerModeTable[2];

static jclass searcherClass;
static jclass cdkExceptionClass;
static jclass inchiExceptionClass;
static jclass tautomerExceptionClass;
static jmethodID getMessageMethod;
static jmethodID getMethod;
static jmethodID indexSizeMethod;
static jmethodID subsearchMethod;
static jmethodID simsearchMethod;
static jfieldID nameField;
static jfieldID lengthField;
static jfieldID idsField;
static jfieldID scoresField;


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


    /* radical modes */
    {
        Oid typoid = LookupExplicitEnumType(spaceid, "radical_mode");

        jclass enumClass = (*env)->FindClass(env, "cz/iocb/sachem/molecule/RadicalMode");
        java_check_exception(__func__);

        const char* const values[] = { "IGNORE", "DEFAULT_AS_STANDARD", "DEFAULT_AS_ANY" };

        for(int i = 0; i < 3; i++)
        {
            radicalModeTable[i].oid = LookupExplicitEnumValue(typoid, values[i]);
            radicalModeTable[i].object = LookupJavaEnumValue(enumClass, values[i], "Lcz/iocb/sachem/molecule/RadicalMode;");
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
            #if PG_VERSION_NUM >= 120000
            desc = CreateTemplateTupleDesc(2);
            #else
            desc = CreateTemplateTupleDesc(2, false);
            #endif

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
    if(unlikely(configQueryPlan == NULL))
    {
        if(unlikely(SPI_connect() != SPI_OK_CONNECT))
            elog(ERROR, "%s: SPI_connect() failed", __func__);

        SPIPlanPtr plan = SPI_prepare("select version, threads from sachem.configuration where index_name = $1", 1, (Oid[]) { VARCHAROID });

        if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        configQueryPlan = plan;

        SPI_finish();
    }


    /* java init */
    searcherClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "cz/iocb/sachem/lucene/Searcher"));
    java_check_exception(__func__);

    cdkExceptionClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "org/openscience/cdk/exception/CDKException"));
    java_check_exception(__func__);

    inchiExceptionClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "cz/iocb/sachem/molecule/InChITools$InChIException"));
    java_check_exception(__func__);

    tautomerExceptionClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "cz/iocb/sachem/molecule/InChITautomerGenerator$InChITautomerException"));
    java_check_exception(__func__);

    getMessageMethod = (*env)->GetMethodID(env, cdkExceptionClass, "getMessage", "()Ljava/lang/String;");
    java_check_exception(__func__);

    getMethod = (*env)->GetStaticMethodID(env, searcherClass, "get", "(Ljava/lang/String;Ljava/lang/String;I)Lcz/iocb/sachem/lucene/Searcher;");
    java_check_exception(__func__);

    indexSizeMethod = (*env)->GetMethodID(env, searcherClass, "indexSize", "()I");
    java_check_exception(__func__);

    subsearchMethod = (*env)->GetMethodID(env, searcherClass, "subsearch", "([BIZLcz/iocb/sachem/molecule/SearchMode;Lcz/iocb/sachem/molecule/ChargeMode;Lcz/iocb/sachem/molecule/IsotopeMode;Lcz/iocb/sachem/molecule/RadicalMode;Lcz/iocb/sachem/molecule/StereoMode;Lcz/iocb/sachem/molecule/AromaticityMode;Lcz/iocb/sachem/molecule/TautomerMode;J)Lcz/iocb/sachem/lucene/SearchResult;");
    java_check_exception(__func__);

    simsearchMethod = (*env)->GetMethodID(env, searcherClass, "simsearch", "([BIZFILcz/iocb/sachem/molecule/AromaticityMode;Lcz/iocb/sachem/molecule/TautomerMode;)Lcz/iocb/sachem/lucene/SearchResult;");
    java_check_exception(__func__);


    jclass resultClass = (*env)->FindClass(env, "cz/iocb/sachem/lucene/SearchResult");
    java_check_exception(__func__);

    nameField = (*env)->GetFieldID(env, resultClass, "name", "Ljava/lang/String;");
    java_check_exception(__func__);

    lengthField = (*env)->GetFieldID(env, resultClass, "length", "I");
    java_check_exception(__func__);

    idsField = (*env)->GetFieldID(env, resultClass, "ids", "[I");
    java_check_exception(__func__);

    scoresField = (*env)->GetFieldID(env, resultClass, "scores", "[F");
    java_check_exception(__func__);


    initialized = true;
}


static jobject lucene_get(VarChar *index)
{
    lucene_search_init();


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    if(unlikely(SPI_execute_plan(configQueryPlan, (Datum[]) { PointerGetDatum(index) }, NULL, true, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 2))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    int32 version = DatumGetInt32(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1));
    int32 threads = DatumGetInt32(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2));

    SPI_finish();

    if(version == 0)
        elog(ERROR, "index has not been synced yet");


    jstring name = NULL;
    jstring folder = NULL;
    jobject instance = NULL;

    PG_TRY();
    {
        name = (*env)->NewStringUTF(env, text_to_cstring(index));
        java_check_exception(__func__);

        folder = (*env)->NewStringUTF(env, get_index_path(text_to_cstring(index), version));
        java_check_exception(__func__);

        instance = (*env)->CallStaticObjectMethod(env, searcherClass, getMethod, name, folder, threads);
        java_check_exception(__func__);

        JavaDeleteRef(name);
    }
    PG_CATCH();
    {
        JavaDeleteRef(name);
        JavaDeleteRef(folder);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return instance;
}


static void lucene_free(jobject lucene)
{
    JavaDeleteRef(lucene);
}


static int32 lucene_index_size(jobject lucene)
{
    int32 size = (*env)->CallIntMethod(env, lucene, indexSizeMethod);
    java_check_exception(__func__);
    return size;
}


static HeapTuple lucene_result_get_item(LuceneResult *result)
{
    while(result->possition < result->length && result->scores[result->possition] == 0)
    {
        char *idx = text_to_cstring(result->index);
        const char *name = (*env)->GetStringUTFChars(env, result->name, NULL);

        if(name == NULL)
            elog(NOTICE, "<unknown>: isomorphism: iteration limit exceeded for target %i in index '%s'", result->ids[result->possition++], idx);
        else if(name[0] == '\0')
            elog(NOTICE, "<unnamed>: isomorphism: iteration limit exceeded for target %i in index '%s'", result->ids[result->possition++], idx);
        else
            elog(NOTICE, "'%s': isomorphism: iteration limit exceeded for target %i in index '%s'", name, result->ids[result->possition++], idx);

        pfree(idx);

        if(name != NULL)
            (*env)->ReleaseStringUTFChars(env, result->name, name);
    }

    if(unlikely(result->possition == result->length))
        return NULL;

    bool isnull[2] = {0, 0};
    Datum values[2] = {Int32GetDatum(result->ids[result->possition]), Float4GetDatum(result->scores[result->possition])};
    HeapTuple tuple = heap_form_tuple(tupdesc, values, isnull);

    result->possition++;

    return tuple;
}


static void lucene_result_free(LuceneResult *result)
{
    if(likely(result != NULL))
    {
        JavaDeleteRef(result->name);
        JavaDeleteIntegerArray(result->idsArray, result->ids, JNI_ABORT);
        JavaDeleteFloatArray(result->scoresArray, result->scores, JNI_ABORT);
    }
}


static LuceneResult *lucene_subsearch(jobject lucene, VarChar *index, VarChar *query, int32 topn, bool sort, Oid search,
        Oid charge, Oid isotope, Oid radical, Oid stereo, Oid aromaticity, Oid tautomers, int64 matchingLimit)
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

        handler = (*env)->CallObjectMethod(env, lucene, subsearchMethod, queryArray, topn, sort,
                ConvertEnumValue(searchModeTable, search),
                ConvertEnumValue(chargeModeTable, charge),
                ConvertEnumValue(isotopeModeTable, isotope),
                ConvertEnumValue(radicalModeTable, radical),
                ConvertEnumValue(stereoModeTable, stereo),
                ConvertEnumValue(aromaticityModeTable,  aromaticity),
                ConvertEnumValue(tautomerModeTable,  tautomers),
                matchingLimit);

        jthrowable exception = (*env)->ExceptionOccurred(env);

        if(exception != NULL && (*env)->IsInstanceOf(env, exception, tautomerExceptionClass) && tautomers == tautomerModeTable[1].oid)
        {
            jstring message = (jstring)(*env)->CallObjectMethod(env, exception, getMessageMethod);
            const char *mstr = message != NULL ? (*env)->GetStringUTFChars(env, message, NULL) : NULL;

            elog(WARNING, "tautomers cannot be generated: %s", mstr != NULL ? mstr : "unknown jvm error");

            if(mstr != NULL)
                (*env)->ReleaseStringUTFChars(env, message, mstr);

            (*env)->ExceptionClear(env);
            JavaDeleteRef(message);
            JavaDeleteRef(exception);

            handler = (*env)->CallObjectMethod(env, lucene, subsearchMethod, queryArray, topn, sort,
                    ConvertEnumValue(searchModeTable, search),
                    ConvertEnumValue(chargeModeTable, charge),
                    ConvertEnumValue(isotopeModeTable, isotope),
                    ConvertEnumValue(radicalModeTable, radical),
                    ConvertEnumValue(stereoModeTable, stereo),
                    ConvertEnumValue(aromaticityModeTable,  aromaticity),
                    tautomerModeTable[0].object,
                    matchingLimit);

            exception = (*env)->ExceptionOccurred(env);
        }

        if(exception != NULL && (*env)->IsInstanceOf(env, exception, inchiExceptionClass) && stereo == stereoModeTable[1].oid)
        {
            jstring message = (jstring)(*env)->CallObjectMethod(env, exception, getMessageMethod);
            const char *mstr = message != NULL ? (*env)->GetStringUTFChars(env, message, NULL) : NULL;

            elog(WARNING, "stereo cannot be determined: %s", mstr != NULL ? mstr : "unknown jvm error");

            if(mstr != NULL)
                (*env)->ReleaseStringUTFChars(env, message, mstr);

            (*env)->ExceptionClear(env);
            JavaDeleteRef(message);
            JavaDeleteRef(exception);

            handler = (*env)->CallObjectMethod(env, lucene, subsearchMethod, queryArray, topn, sort,
                    ConvertEnumValue(searchModeTable, search),
                    ConvertEnumValue(chargeModeTable, charge),
                    ConvertEnumValue(isotopeModeTable, isotope),
                    ConvertEnumValue(radicalModeTable, radical),
                    stereoModeTable[0].object,
                    ConvertEnumValue(aromaticityModeTable,  aromaticity),
                    tautomerModeTable[0].object,
                    matchingLimit);
        }

        java_check_exception(__func__);

        JavaDeleteRef(queryArray);

        result = palloc0(sizeof(LuceneResult));

        result->index = index;
        result->name = (jstring) (*env)->GetObjectField(env, handler, nameField);

        result->idsArray = (*env)->GetObjectField(env, handler, idsField);
        result->ids = (*env)->GetIntArrayElements(env, result->idsArray, NULL);

        result->scoresArray = (*env)->GetObjectField(env, handler, scoresField);
        result->scores = (*env)->GetFloatArrayElements(env, result->scoresArray, NULL);

        result->length = (*env)->GetIntField(env, handler, lengthField);

        if(result->ids == NULL || result->scores == NULL)
            elog(ERROR, "out of memmory");;
    }
    PG_CATCH();
    {
        JavaDeleteRef(queryArray);
        JavaDeleteRef(handler);
        lucene_result_free(result);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return result;
}


static LuceneResult *lucene_simsearch(jobject lucene, VarChar *index, VarChar *query, int32 topn, bool sort,
        float4 threshold, int32 radius, Oid aromaticity, Oid tautomers)
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

        handler = (*env)->CallObjectMethod(env, lucene, simsearchMethod, queryArray, topn, sort, threshold, radius,
                ConvertEnumValue(aromaticityModeTable, aromaticity),
                ConvertEnumValue(tautomerModeTable, tautomers));

        jthrowable exception = (*env)->ExceptionOccurred(env);

        if(exception != NULL && (*env)->IsInstanceOf(env, exception, tautomerExceptionClass) && tautomers == tautomerModeTable[1].oid)
        {
            jstring message = (jstring)(*env)->CallObjectMethod(env, exception, getMessageMethod);
            const char *mstr = message != NULL ? (*env)->GetStringUTFChars(env, message, NULL) : NULL;

            elog(WARNING, "tautomers cannot be generated: %s", mstr != NULL ? mstr : "unknown jvm error");

            if(mstr != NULL)
                (*env)->ReleaseStringUTFChars(env, message, mstr);

            (*env)->ExceptionClear(env);
            JavaDeleteRef(message);
            JavaDeleteRef(exception);

            handler = (*env)->CallObjectMethod(env, lucene, simsearchMethod, queryArray, topn, sort, threshold, radius,
                    ConvertEnumValue(aromaticityModeTable, aromaticity),
                    tautomerModeTable[0].object);
        }

        java_check_exception(__func__);

        JavaDeleteRef(queryArray);

        result = palloc0(sizeof(LuceneResult));

        result->index = index;
        result->name = (jstring) (*env)->GetObjectField(env, handler, nameField);

        result->idsArray = (*env)->GetObjectField(env, handler, idsField);
        result->ids = (*env)->GetIntArrayElements(env, result->idsArray, NULL);

        result->scoresArray = (*env)->GetObjectField(env, handler, scoresField);
        result->scores = (*env)->GetFloatArrayElements(env, result->scoresArray, NULL);

        result->length = (*env)->GetArrayLength(env, result->idsArray);

        if(result->ids == NULL || result->scores == NULL)
            elog(ERROR, "out of memmory");;
    }
    PG_CATCH();
    {
        JavaDeleteRef(queryArray);
        JavaDeleteRef(handler);
        lucene_result_free(result);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return result;
}


PG_FUNCTION_INFO_V1(index_size);
Datum index_size(PG_FUNCTION_ARGS)
{
    VarChar *index = PG_GETARG_VARCHAR_P(0);

    jobject lucene = lucene_get(index);
    int32 size;

    PG_TRY();
    {
        size = lucene_index_size(lucene);

        lucene_free(lucene);
    }
    PG_CATCH();
    {
        lucene_free(lucene);
        PG_RE_THROW();
    }
    PG_END_TRY();

    return size;
}


PG_FUNCTION_INFO_V1(substructure_search);
Datum substructure_search(PG_FUNCTION_ARGS)
{
    if(unlikely(SRF_IS_FIRSTCALL()))
    {
        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();

        VarChar *index = PG_GETARG_VARCHAR_P(0);
        VarChar *query = PG_GETARG_VARCHAR_P(1);
        Oid search = PG_GETARG_OID(2);
        Oid charge = PG_GETARG_OID(3);
        Oid isotope = PG_GETARG_OID(4);
        Oid radical = PG_GETARG_OID(5);
        Oid stereo = PG_GETARG_OID(6);
        Oid aromaticity = PG_GETARG_OID(7);
        Oid tautomers = PG_GETARG_OID(8);
        int32 topn = PG_GETARG_INT32(9);
        bool sort = PG_GETARG_BOOL(10);
        int64 matchingLimit = PG_GETARG_INT64(11);

        jobject lucene = lucene_get(index);

        PG_TRY();
        {
            PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);
            VarChar *index = PG_GETARG_VARCHAR_P(0);
            funcctx->user_fctx = lucene_subsearch(lucene, index, query, topn, sort, search, charge, isotope, radical, stereo, aromaticity, tautomers, matchingLimit);
            PG_MEMCONTEXT_END();

            lucene_free(lucene);
        }
        PG_CATCH();
        {
            lucene_free(lucene);
            PG_RE_THROW();
        }
        PG_END_TRY();
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

    if(likely(HeapTupleIsValid(item)))
        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(item));

    lucene_result_free(result);
    SRF_RETURN_DONE(funcctx);
}


PG_FUNCTION_INFO_V1(similarity_search);
Datum similarity_search(PG_FUNCTION_ARGS)
{
    if(unlikely(SRF_IS_FIRSTCALL()))
    {
        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();

        VarChar *index = PG_GETARG_VARCHAR_P(0);
        VarChar *query = PG_GETARG_VARCHAR_P(1);
        float4 threshold = PG_GETARG_FLOAT4(2);
        int32 radius = PG_GETARG_INT32(3);
        Oid aromaticity = PG_GETARG_OID(4);
        Oid tautomers = PG_GETARG_OID(5);
        int32 topn = PG_GETARG_INT32(6);
        bool sort = PG_GETARG_BOOL(7);

        jobject lucene = lucene_get(index);

        PG_TRY();
        {
            PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);
            VarChar *index = PG_GETARG_VARCHAR_P(0);
            funcctx->user_fctx = lucene_simsearch(lucene, index, query, topn, sort, threshold, radius, aromaticity, tautomers);
            PG_MEMCONTEXT_END();

            lucene_free(lucene);
        }
        PG_CATCH();
        {
            lucene_free(lucene);
            PG_RE_THROW();
        }
        PG_END_TRY();
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

    if(likely(HeapTupleIsValid(item)))
        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(item));

    lucene_result_free(result);
    SRF_RETURN_DONE(funcctx);
}
