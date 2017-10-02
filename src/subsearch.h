#ifndef SUBSEARCH_H_
#define SUBSEARCH_H_

#define FETCH_SIZE           1000
#define FP_SIZE               788
#define COUNTS_SIZE            10

#define MOLECULES_TABLE          "molecules"
#define MOLECULE_COUNTS_TABLE    "molecule_counts"
#define FINGERPRINT_INDEX_TABLE  "fingerprint_orchem_index"

#define COUNT_INDEX_FILE        "molecule_counts.idx"
#define FINGERPRINT_INDEX_FILE  "fingerprint_orchem_index.idx"

void subsearch_module_init(void);
void subsearch_module_finish(void);

#endif /* SUBSEARCH_H_ */
