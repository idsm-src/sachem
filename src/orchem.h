#ifndef ORCHEM_H_
#define ORCHEM_H_

#define FETCH_SIZE           1000
#define FP_SIZE               788
#define COUNTS_SIZE            10

#define MOLECULES_TABLE          "molecules"
#define MOLECULE_COUNTS_TABLE    "molecule_counts"
#define FINGERPRINT_INDEX_TABLE  "fingerprint_orchem_index"

void orchem_module_init(void);
void orchem_module_finish(void);

#endif /* ORCHEM_H_ */
