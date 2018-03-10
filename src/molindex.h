#ifndef MOLINDEX_H_
#define MOLINDEX_H_

#include <stdbool.h>


#define MOLECULE_INDEX_PREFIX         "sachem_molecules"
#define MOLECULE_INDEX_SUFFIX         ".idx"


void sachem_generate_molecule_index(int indexNumber, bool useSeqId);

#endif /* MOLINDEX_H_ */
