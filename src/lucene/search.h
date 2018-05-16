#ifndef LUCENE_SEARCH_H_
#define LUCENE_SEARCH_H_

#include "lucene.h"


extern Lucene lucene;


void lucene_search_init(void);
int lucene_search_update_snapshot(void);

#endif /* LUCENE_SEARCH_H_ */
