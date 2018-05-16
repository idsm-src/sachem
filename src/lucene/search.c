#include <postgres.h>
#include <executor/spi.h>
#include "common.h"
#include "search.h"
#include "sachem.h"
#include "lucene.h"
#include "java/parse.h"


static bool initialized = false;
static bool javaInitialized = false;
static bool luceneInitialised = false;
static int indexId = -1;
static SPIPlanPtr snapshotQueryPlan;
Lucene lucene;


void lucene_search_init(void)
{
    if(unlikely(initialized == false))
    {
        if(unlikely(javaInitialized == false))
        {
            java_parse_init();
            javaInitialized = true;
        }


        /* prepare lucene */
        if(unlikely(luceneInitialised == false))
        {
            lucene_init(&lucene);
            luceneInitialised = true;
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


        initialized = true;
    }
}


int lucene_search_update_snapshot(void)
{
    if(unlikely(SPI_connect() != SPI_OK_CONNECT))
        elog(ERROR, "%s: SPI_connect() failed", __func__);

    if(unlikely(SPI_execute_plan(snapshotQueryPlan, NULL, NULL, true, 0) != SPI_OK_SELECT))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    if(unlikely(SPI_processed != 1 || SPI_tuptable == NULL || SPI_tuptable->tupdesc->natts != 1))
        elog(ERROR, "%s: SPI_execute_plan() failed", __func__);

    char isNullFlag;
    int32_t dbIndexNumber = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isNullFlag);

    if(unlikely(SPI_result == SPI_ERROR_NOATTRIBUTE || isNullFlag))
        elog(ERROR, "%s: SPI_getbinval() failed", __func__);

    SPI_finish();


    if(unlikely(dbIndexNumber != indexId))
    {
        char *path = get_index_path(LUCENE_INDEX_PREFIX, LUCENE_INDEX_SUFFIX, dbIndexNumber);
        lucene_set_folder(&lucene, path);
        indexId = dbIndexNumber;
    }

    return indexId;
}
