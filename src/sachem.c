#include <postgres.h>
#include <fmgr.h>
#include "isomorphism.h"
#include "java.h"
#include "simsearch.h"
#include "subsearch.h"
#include "subsearch-gin.h"


PG_MODULE_MAGIC;


void __attribute__ ((constructor)) sachemInit(void)
{
    isomorphism_module_init();
    java_module_init();
    simsearch_module_init();
    subsearch_module_init();
    subsearch_gin_module_init();
}


void __attribute__ ((destructor)) sachemFinish(void)
{
    isomorphism_module_finish();
    java_module_finish();
    simsearch_module_finish();
    subsearch_module_finish();
    subsearch_gin_module_finish();
}
