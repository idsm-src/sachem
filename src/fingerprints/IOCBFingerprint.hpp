#ifndef IOCB_FINGERPRINT_HPP__
#define IOCB_FINGERPRINT_HPP__

#include <map>
#include <set>


typedef struct Molecule Molecule;
typedef std::map<uint32_t, std::set<uint32_t>> BitInfo;


std::set<uint32_t> iocb_fingerprint_get(const Molecule *molecule, int graphSize, int ringSize, int maxFeatLogCount,
        bool forQuery = false, BitInfo *info = nullptr);

#endif /* IOCB_FINGERPRINT_HPP__ */
