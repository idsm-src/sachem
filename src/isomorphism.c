#include <postgres.h>
#include "isomorphism.h"


#if USE_VF2_TIMEOUT
TimeoutId vf2TimeoutId = -1;
volatile bool vf2Timeouted;


static void vf2_timeout_handler()
{
    elog(WARNING, "isomorphism: VF2 timeout expired");
    vf2Timeouted = true;
}
#endif


void isomorphism_module_init(void)
{
#if USE_VF2_TIMEOUT
    PG_TRY();
    {
        vf2TimeoutId = RegisterTimeout(USER_TIMEOUT, vf2_timeout_handler);
    }
    PG_CATCH();
    {
        elog(WARNING, "%s: RegisterTimeout() failed", __func__);
    }
    PG_END_TRY();
#endif
}


void isomorphism_module_finish(void)
{
#if USE_VF2_TIMEOUT
    vf2TimeoutId = -1;
#endif
}
