#ifndef CRNG_FINGERPRINT_HPP__
#define CRNG_FINGERPRINT_HPP__

#include <map>
#include <set>


typedef struct Molecule Molecule;
typedef std::map<uint32_t, std::set<uint32_t>> BitInfo;


std::map<uint32_t, int> crng_fingerprint_get(const Molecule *molecule, BitInfo *info = nullptr);

#endif /* CRNG_FINGERPRINT_HPP__ */
