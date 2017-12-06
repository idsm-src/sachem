#include <postgres.h>
#include "isomorphism.h"


#if USE_VF2_TIMEOUT
TimeoutId vf2TimeoutId = MAX_TIMEOUTS;
volatile bool vf2Timeouted;


void vf2_timeout_handler()
{
    elog(WARNING, "isomorphism: VF2 timeout expired");
    vf2Timeouted = true;
}
#endif
