
#include "CRNGFingerprint.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <memory>
#include <boost/thread/mutex.hpp>

namespace RDKit
{
namespace CRNGFingerprints
{

//some most common rings from ChEBI
static const char* patterns[] = {
	"C1=CCCCC=CCCC1",
	"C1=CCCCCCCCCCC1",
	"C1CCCCCOCCCCC1",
	"C1CCCCCOCCCOCCOCCCCC1",
	"C1=CCCOCC1",
	"C1=CC[NH2+]CC1",
	"C1CNPOC1",
	"C1COCCCOCCCOC1",
	"C1COCC[NH2+]1",
	"c1nCCCOPOCCNCCCCCN[Co-3][nH+]1",
	"c1nCCCS1",
	"C1=NCNCN1",
	"C1NCNN1",
	"O1[W]O[W]O[W]1",
	"c1ccc1",
	"C1=CCCCCCCC1",
	"C1=CCCCCCC=CCC1",
	"C1ccCCN=1",
	"C1cc-ncCN=1",
	"c1Cc[nH]cCc[nH]cCc[nH]cCc[nH]1",
	"C1=CCOCCCCC1",
	"c1ccOccccCCCNCCccccOcccCCNCCC1",
	"C1CNCCNCCNCCNCCNCCSSC1",
	"c1-cn[nH]c1",
	"c1cOccCC1",
	"C1NCCC[NH+]=1",
	"C1NCC[NH+]=1",
	"C1NCSCN1",
	"O1[W]O[W]O[W]O[W]1",
	"c1cCCCCCC1",
	"C1=CCCC=CCCCC=CCC1",
	"C1=CC=CCCCOCCCNCCOCCCCC=CCCCCCCC=C1",
	"C1CCccN=1",
	"C1C=COCC=1",
	"C1CCOCOCCOCOCCOCOCCOCCCCOCOCCOCOCCOCOCCOC1",
	"C1=CN=CC=CNCCN=CC=CN=C1",
	"c1-cncn1",
	"C1CNCN=1",
	"C1=CN[Co-3][NH+]=C1",
	"C1C[NH2+][Co-3]N1",
	"c1cn[nH]cn1",
	"c1cSNC1",
	"c1nnn[nH]1",
	"S1[Fe]S[Fe+]1",
	"c1cCC[nH+]c-1",
	"c1ccOcccCNCCNCC1",
	"c1cCOCCOCCCOCCCC1",
	"c1cNCCNC1",
	"c1-cnc[n-]c1",
	"C1=CNC=NC1",
	"c1cncoc1",
	"C1=C[NH2+][Co-2]nc1",
	"c1cnnnc1",
	"c1cSccC1",
	"c1cSNCN1",
	"c1=NC=CCn1",
	"C1=CCccc1",
	"c1-ccCCCc1",
	"C1CCCCCCCCC1",
	"C1=CCCCCCCCCC=CCC1",
	"c1cCNCCCC1",
	"c1cc[nH]nc1",
	"C1CNCCNCCNCCN1",
	"C1CN=CN=1",
	"c1c[NH2+]CC1",
	"c1c[nH+]c[nH]1",
	"c1cSCC1",
	"c1c[s+]ccn1",
	"[BH2-]1ncC=C[NH2+]1",
	"c1-cccc1",
	"C1=CCCCCCCCCC1",
	"c1cCCNCCCC1",
	"c1-ccCOCCCOCc1",
	"c1=CC=NC=Ccn[Mg]n1",
	"C1CNCCNCCNCCNCCOC1",
	"C1C[NH2+][Gd-][NH2+]1",
	"c1c[nH]ccn1",
	"c1csc[nH]1",
	"c1nCCS1",
	"c1nn[nH+]n1",
	"O1[SiH2]O[SiH2]O[SiH2]O[SiH2]1",
	"C1=CCCCCCC1",
	"c1cC=CC[NH+]=1",
	"C1=CCCCOCCCNCCNCCNCCC1",
	"c1c[n-]cn1",
	"c1cnCOCn1",
	"c1cOCCOCCCNCCCNC1",
	"c1cSCCC1",
	"c1ncsn1",
	"c1cCCnc-1",
	"c1-ccNCc1",
	"C1CCNCCNCCNCCNCCNCCNCCNC1",
	"c1cNccNC1",
	"c1cOccOC1",
	"C1CSCCN1",
	"c1cCccc-1",
	"c1c[nH+]ccn1",
	"c1CCcn1",
	"C1CCCOCCC1",
	"c1ccNC1",
	"c1cCNCCC1",
	"C1=CNccC1",
	"C1CNCOC1",
	"C1CNOC1",
	"c1cOcccCCcccOccccCC1",
	"c1cOccO1",
	"C1=CSCN1",
	"C1CScnN=1",
	"C1C[Zr]1",
	"C1=NccCN1",
	"C1=CCCCCC=CCCC1",
	"C1CCCCOCCC1",
	"c1cC[CH-]C1",
	"c1cCnc-1",
	"C1=CN=CC1",
	"C1C[NH2+]CCN1",
	"C1COOC1",
	"C1=NCNN1",
	"C1NO[Fe-3][O+]=1",
	"c1=CCCc1",
	"C1CCCOCCCCOCC1",
	"c1c[o+]ccn1",
	"c1csnn1",
	"C1NCNCN1",
	"c1c-ccc-1",
	"C1C=CCNcccccCCOC=CCCCCCCCCC=1",
	"c1cNccN1",
	"c1cSCN1",
	"C1=CS[Mo]S1",
	"c1cCOCCOCCCOCCC1",
	"C1C[NH2+][Co-2]N1",
	"C1=C[NH2+][Co-3][NH+]=C1",
	"C1CNNC1",
	"c1nn[nH]n1",
	"c1-ccOCc1",
	"c1cNCN1",
	"c1cOCCCOCCNCCNC1",
	"c1-cnccc1",
	"C1=CN=CNC1",
	"C1=C[NH2+][Co-2][NH+]=C1",
	"c1c[nH]c[nH]1",
	"c1cOcccCCccccOcccCC1",
	"c1coc[nH+]1",
	"c1cSCCN1",
	"C1CCOCOCCOCOCCOCOCCOCOCCOCOCCOCOCCOCOCCOC1",
	"C1=C[NH+]=CC1",
	"c1cOCCCO1",
	"c1cOCN1",
	"C1=CCcc1",
	"C1C=CCC=1",
	"c1CCCCCn1",
	"c1c[nH][nH]c1",
	"C1CSSCN1",
	"c1ccCC1",
	"c1cccccc1",
	"C1CCCNCCC1",
	"c1=CccOc1",
	"c1c=NCC=1",
	"c1cNccCC1",
	"c1-cnccn1",
	"c1CCCCNCCOCn1",
	"C1C=CNCC=1",
	"C1ccNCCN=1",
	"c1cOCCOCCCNCCNC1",
	"C1C[U]1",
	"c1cCCC[NH+]=1",
	"c1-coccn1",
	"C1COPO1",
	"c1-ccCCCCc1",
	"c1cOCCCOCCNCCCNC1",
	"C1CCCCN=1",
	"C1CCNCN=1",
	"c1cCNCCNCCccc-1",
	"C1CSSC1",
	"c1ccCNC1",
	"C1CNCCNC1",
	"c1-cnc[nH]c1",
	"C1=CNNC1",
	"S1[Fe]S[Fe]1",
	"C1CCOOC1",
	"C1COCCO1",
	"c1cCCccC1",
	"c1cCCCCC1",
	"C1CNCCNCCNCCNCCNCCOC1",
	"c1cnc[nH+]c1",
	"c1CCCCn1",
	"C1C=CCCOCCCCCCCCCC=1",
	"c1cnsn1",
	"C1C=CCCCOCCCCC=CCCC=1",
	"C1=CCOCCC1",
	"C1=CN=CCC1",
	"C1CCON=1",
	"c1cccc1",
	"C1CCCCCCC1",
	"c1cnCCCCNCCOC1",
	"C1COC1",
	"c1=CC=[NH+][Fe-2]n1",
	"C1=CN[Co-2][NH+]=C1",
	"c1-cncnc1",
	"C1=CCCCCC=CCC1",
	"C1CCCCCCOCCCCCC1",
	"C1CNccN=1",
	"C1=CNCNC1",
	"C1Ccc[NH+]=1",
	"C1CN1",
	"C1=C[NH2+][Fe-2]nc1",
	"C1=COccC1",
	"C1=NCCCN1",
	"C1CccN=1",
	"C1=CNCCNCCCOcccc1",
	"C1=CCNcc1",
	"C1=COC=CC1",
	"c1cnccc[nH]cccnccc[nH]c1",
	"c1cNCNC1",
	"C1=CCCCCCCCC1",
	"c1cCCNCC1",
	"c1cnsc1",
	"c1cOccC1",
	"c1-ccCc1",
	"C1NCNCNCN1",
	"c1ccOccccCCCNCCNC1",
	"c1csc[nH+]1",
	"c1cOCCN1",
	"C1COPOC1",
	"C1CCSC1",
	"c1cnncn1",
	"c1cnon1",
	"c1ccOccccCCNCCNCC1",
	"C1=COCC1",
	"C1=CNCC1",
	"c1-coccc1",
	"C1=CNCCC1",
	"C1=NCCO1",
	"c1cC[NH2+]CC1",
	"c1ccCCC1",
	"c1c[nH+]cn1",
	"C1C[Fe]1",
	"c1cSccN1",
	"c1ncon1",
	"c1cOCCCOCCNCCCCNC1",
	"C1CCNN=1",
	"c1c[nH]c[nH]c1",
	"C1=CNC=CC1",
	"C1CCCOCC1",
	"c1cNCCN1",
	"C1=CCN=C1",
	"c1cOCCOCCCNCCCCNC1",
	"C1=CCCCCC1",
	"C1CCCNCC1",
	"c1-ccCCc1",
	"c1-ccccc1",
	"c1=CC=[NH+][MgH2-2]n1",
	"C1=C[NH2+][MgH2-2]nc1",
	"C1COCN1",
	"C1=CC[NH+]=C1",
	"C1=NCCN1",
	"C1COCOC1",
	"C1=NCCS1",
	"c1cCOC1",
	"C1CCCN=1",
	"C1=CNCSC1",
	"C1CC[NH2+]C1",
	"c1cCOCC1",
	"C1=CCNCC1",
	"C1=CCCcc1",
	"c1nc[nH]n1",
	"C1=CCOcc1",
	"C1CCC[NH+]=1",
	"C1COCO1",
	"c1ccnnc1",
	"C1=CCNC1",
	"C1CNCNC1",
	"C1CC[NH2+]CC1",
	"c1nnco1",
	"c1ncncn1",
	"C1=CCccC1",
	"c1cc[o+]cc1",
	"c1cn[nH]c1",
	"c1nncs1",
	"C1=COCCC1",
	"C1CCCCCC1",
	"c1cOCCO1",
	"c1cc[nH]cc1",
	"c1cocn1",
	"c1CCCn1",
	"c1nnnn1",
	"c1CNCCn1",
	"c1cc[nH+]cc1",
	"c1cCNC1",
	"C1CCNCCNC1",
	"c1cCCC1",
	"c1cCccC1",
	"C1CNCN1",
	"C1C=CCCC=1",
	"C1CSCN1",
	"C1=CCOC1",
	"c1COCCCNCCCCn1",
	"C1=CCC=CC1",
	"c1cnCCCCNCCCOC1",
	"c1ncnn1",
	"c1cNCC1",
	"C1CCC1",
	"C1CNCCCOC1",
	"c1cnccn1",
	"c1cNCCC1",
	"c1c[nH]cn1",
	"C1=CCCC1",
	"c1cOCCC1",
	"c1-ccCOCCCNCc1",
	"c1cnnc1",
	"C1CO1",
	"c1cnoc1",
	"c1cOCCCNCC1",
	"c1cCCCC1",
	"C1=CCOCC1",
	"c1cSNCCCO1",
	"C1COCCN1",
	"c1cscn1",
	"c1cCNCC1",
	"c1cnnn1",
	"c1ccocc1",
	"c1ccsc1",
	"c1ccoc1",
	"c1cnc[nH]c1",
	"c1ccnc1",
	"c1cOCO1",
	"c1cc[nH]c1",
	"c1cOCC1",
	"C1CNCCN1",
	"c1cOCCCCCOCCCNC1",
	"C1CC1",
	"c1cOCCNCCCCNC1",
	"C1CCNC1",
	"C1CNC1",
	"c1cOCCCNC1",
	"C1CCNCC1",
	"c1cncn1",
	"C1=CCCCC1",
	"c1cncnc1",
	"C1CCCC1",
	"C1CCOC1",
	"c1ccncc1",
	"C1CCCCC1",
	"C1CCOCC1",
	"c1ccccc1",
	nullptr
};

static std::vector<std::unique_ptr<ROMol>> smarts;
static bool init = false;
static boost::mutex mtx;

SparseIntVect<boost::uint32_t> *getFingerprint (const ROMol &mol,
                                                BitInfo *info)
{
	if (!init) {
	    mtx.lock();
	    if (!init) {
	        size_t cnt = 0;
	        //count the patterns
	        for (const char**pp = patterns; *pp; ++pp, ++cnt);
	        smarts.resize (cnt);
	        for (size_t i = 0; i < cnt; ++i)
	            smarts[i] = std::unique_ptr<ROMol>(SmartsToMol (patterns[i]));
	        init = true;
	    }
		mtx.unlock();
	}

	std::map<uint32_t, int> fp;
	std::vector<MatchVectType> matches;
	for (size_t i = 0; i < smarts.size(); ++i) {
		int nmatches = SubstructMatch
		               (mol, *smarts[i], matches,
		                true, false, false, false,
		                256); //lower could probably suffice
		if (nmatches) {
			fp[i] = nmatches;
			if (info)
				for (auto&mv : matches)
					for (auto&m : mv)
						(*info) [fp[i]]
						.insert (m.second);
		}
	}

	SparseIntVect<boost::uint32_t> *res
	    = new SparseIntVect<boost::uint32_t> (smarts.size() + 1);
	for (auto&i : fp) res->setVal (i.first, i.second);
	return res;
}

}
}
