#ifndef RC_FINGERPRINT_HPP__
#define RC_FINGERPRINT_HPP__

#include <map>
#include <set>


typedef struct Molecule Molecule;
typedef std::map<uint32_t, std::set<uint32_t>> BitInfo;


std::map<uint32_t, int> rc_fingerprint_get(const Molecule *molecule, int minRadius, int maxRadius,
        BitInfo *info = nullptr);

#endif /* RC_FINGERPRINT_HPP__ */
