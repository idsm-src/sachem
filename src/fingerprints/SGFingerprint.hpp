#ifndef SG_FINGERPRINT_HPP__
#define SG_FINGERPRINT_HPP__

#include <map>
#include <set>


typedef struct Molecule Molecule;
typedef std::map<uint32_t, std::set<uint32_t>> BitInfo;


std::map<uint32_t, int> sg_fingerprint_get(const Molecule *molecule, int lowerLen, int upperLen, bool forQuery = false,
        BitInfo *info = nullptr);

#endif /* SG_FINGERPRINT_HPP__ */
