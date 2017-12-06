#include <postgres.h>
#include <fmgr.h>
#include "isomorphism.h"
#include "java.h"
#include "subsearch-lucy.h"


PG_MODULE_MAGIC;


void __attribute__ ((constructor)) sachemInit(void)
{
    java_module_init();
    subsearch_lucy_module_init();
}


void __attribute__ ((destructor)) sachemFinish(void)
{
    java_module_finish();
    subsearch_lucy_module_finish();
}
