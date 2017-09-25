#include <postgres.h>
#include <fmgr.h>
#include "java.h"
#include "orchem.h"
#include "simsearch.h"


PG_MODULE_MAGIC;


void __attribute__ ((constructor)) pgchemInit(void)
{
    java_module_init();
    orchem_module_init();
    simsearch_module_init();
}


void __attribute__ ((destructor)) pgchemFinish(void)
{
    java_module_finish();
    orchem_module_finish();
    simsearch_module_finish();
}
