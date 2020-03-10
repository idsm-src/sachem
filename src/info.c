#include <postgres.h>
#include <access/htup_details.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <utils/memutils.h>
#include <funcapi.h>
#include "java.h"
#include "sachem.h"


static bool initialized = false;
static TupleDesc tupdesc = NULL;
static SPIPlanPtr configQueryPlan = NULL;

static jclass indexInfoClass;
static jmethodID getSegmentInfosMethod;
static jfieldID nameField;
static jfieldID moleculesField;
static jfieldID deletesField;
static jfieldID sizeField;


static void lucene_info_init(void)
{
    if(likely(initialized))
        return;

    java_init();


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
            desc = CreateTemplateTupleDesc(4);
            #else
            desc = CreateTemplateTupleDesc(4, false);
            #endif

            TupleDescInitEntry(desc, (AttrNumber) 1, "name", VARCHAROID, -1, 0);
            TupleDescInitEntry(desc, (AttrNumber) 2, "molecules", INT4OID, -1, 0);
            TupleDescInitEntry(desc, (AttrNumber) 3, "deletes", INT4OID, -1, 0);
            TupleDescInitEntry(desc, (AttrNumber) 4, "size", INT8OID, -1, 0);
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

        SPIPlanPtr plan = SPI_prepare("select version from sachem.configuration where index_name = $1", 1, (Oid[]) { VARCHAROID });

        if(unlikely(SPI_keepplan(plan) == SPI_ERROR_ARGUMENT))
            elog(ERROR, "%s: SPI_keepplan() failed", __func__);

        configQueryPlan = plan;

        SPI_finish();
    }


    /* java init */
    indexInfoClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "cz/iocb/sachem/lucene/IndexInfo"));
    java_check_exception(__func__);

    getSegmentInfosMethod = (*env)->GetStaticMethodID(env, indexInfoClass, "getSegmentInfos", "(Ljava/lang/String;)[Lcz/iocb/sachem/lucene/IndexInfo$SegmentInfo;");
    java_check_exception(__func__);

    jclass segmentInfoClass = (*env)->FindClass(env, "cz/iocb/sachem/lucene/IndexInfo$SegmentInfo");
    java_check_exception(__func__);

    nameField = (*env)->GetFieldID(env, segmentInfoClass, "name", "Ljava/lang/String;");
    java_check_exception(__func__);

    moleculesField = (*env)->GetFieldID(env, segmentInfoClass, "molecules", "I");
    java_check_exception(__func__);

    deletesField = (*env)->GetFieldID(env, segmentInfoClass, "deletes", "I");
    java_check_exception(__func__);

    sizeField = (*env)->GetFieldID(env, segmentInfoClass, "size", "J");
    java_check_exception(__func__);


    initialized = true;
}


static HeapTuple *lucene_get_segment_infos(VarChar *index)
{
    lucene_info_init();


    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    if(unlikely(SPI_execute_plan(configQueryPlan, (Datum[]) { PointerGetDatum(index) }, NULL, true, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    int32 version = DatumGetInt32(SPI_get_value(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1));

    SPI_finish();


    jstring folder = NULL;
    jarray infos = NULL;
    jobject item = NULL;
    jobject name = NULL;

    HeapTuple *result;

    PG_TRY();
    {
        folder = (*env)->NewStringUTF(env, get_index_path(text_to_cstring(index), version));
        java_check_exception(__func__);

        infos = (*env)->CallStaticObjectMethod(env, indexInfoClass, getSegmentInfosMethod, folder);
        java_check_exception(__func__);

        jsize count = (*env)->GetArrayLength(env, infos);


        result = palloc((count + 1) * sizeof(HeapTuple));
        result[count] = NULL;

        for(int i = 0; i < count; i++)
        {
            item = (*env)->GetObjectArrayElement(env, infos, i);

            name = (*env)->GetObjectField(env, item, nameField);
            const char *cname = (*env)->GetStringUTFChars(env, name, NULL);
            java_check_exception(__func__);

            JavaDeleteRef(name);

            jint molecules = (*env)->GetIntField(env, item, moleculesField);
            jint deletes = (*env)->GetIntField(env, item, deletesField);
            jlong size = (*env)->GetLongField(env, item, sizeField);

            JavaDeleteRef(item);

            bool isnull[4] = {0, 0, 0, 0};
            Datum values[4] = {PointerGetDatum(cstring_to_text(cname)), Int32GetDatum(molecules), Int32GetDatum(deletes), Int64GetDatum(size)};
            result[i] = heap_form_tuple(tupdesc, values, isnull);
        }

        JavaDeleteRef(folder);
        JavaDeleteRef(infos);
    }
    PG_CATCH();
    {
        JavaDeleteRef(folder);
        JavaDeleteRef(infos);
        JavaDeleteRef(item);
        JavaDeleteRef(name);

        PG_RE_THROW();
    }
    PG_END_TRY();

    return result;
}


PG_FUNCTION_INFO_V1(segments);
Datum segments(PG_FUNCTION_ARGS)
{
    if(unlikely(SRF_IS_FIRSTCALL()))
    {
        FuncCallContext *funcctx = SRF_FIRSTCALL_INIT();
        VarChar *index = PG_GETARG_VARCHAR_P(0);

        PG_MEMCONTEXT_BEGIN(funcctx->multi_call_memory_ctx);
        funcctx->user_fctx = lucene_get_segment_infos(index);
        PG_MEMCONTEXT_END();
    }


    FuncCallContext *funcctx = SRF_PERCALL_SETUP();
    HeapTuple *tuples = (HeapTuple *) funcctx->user_fctx;
    HeapTuple tuple = tuples[funcctx->call_cntr];

    if(likely(HeapTupleIsValid(tuple)))
        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(tuple));

    SRF_RETURN_DONE(funcctx);
}
