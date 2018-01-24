
#ifndef __FPSEARCH_IOCBFPS_H__
#define __FPSEARCH_IOCBFPS_H__

#include <DataStructs/SparseIntVect.h>
#include <boost/cstdint.hpp>
#include <map>
#include <set>

namespace RDKit
{
class ROMol;
namespace IOCBFingerprints
{
const std::string version = "0.0.1";

typedef std::map<boost::uint32_t, std::set<boost::uint32_t> > BitInfo;

//! returns the IOCB fingerprint for a molecule
/*!
  \param mol:    the molecule to be fingerprinted

  \return a pointer to the fingerprint. The client is
  responsible for calling delete on this.

*/
SparseIntVect<boost::uint32_t> *getFingerprint (const ROMol &mol,
                                                int graphSize,
                                                int ringSize,
                                                int maxFeatLogCount,
                                                bool forQuery = false,
                                                BitInfo *info = nullptr);
}
}

#endif
