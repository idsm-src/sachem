#ifndef ENVIROMENT_H_
#define ENVIROMENT_H_

#include <postgres.h>
#include <miscadmin.h>


/*
 * Tells whether a pending interrupt is going to end the current query. InterruptPending alone does not say so, because
 * PostgreSQL also sets it for interrupts that are merely processed, e.g. the flush of statistics 10 s after a short
 * statement. The flags can be read in any thread.
 */
#if PG_VERSION_NUM >= 170000
#define QUERY_CANCEL_PENDING() (unlikely(InterruptPending) && (QueryCancelPending || ProcDiePending || ClientConnectionLost || TransactionTimeoutPending))
#else
#define QUERY_CANCEL_PENDING() (unlikely(InterruptPending) && (QueryCancelPending || ProcDiePending || ClientConnectionLost))
#endif

#endif /* ENVIROMENT_H_ */
