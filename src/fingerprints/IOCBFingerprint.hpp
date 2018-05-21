#ifndef IOCB_FINGERPRINT_HPP__
#define IOCB_FINGERPRINT_HPP__

#include <map>
#include <set>


typedef struct Molecule Molecule;
typedef std::map<uint32_t, std::set<uint32_t>> BitInfo;


std::set<uint32_t> iocb_substructure_fingerprint_get(const Molecule *molecule, int graphSize, int maxFeatLogCount,
        bool forQuery = false, BitInfo *info = nullptr);
std::set<uint32_t> iocb_similarity_fingerprint_get(const Molecule *molecule, int circSize, int maxFeatLogCount,
        BitInfo *info = nullptr);

#endif /* IOCB_FINGERPRINT_HPP__ */
