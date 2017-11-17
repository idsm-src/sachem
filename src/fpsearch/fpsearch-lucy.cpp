#include "fpsearch-lucy.h"

#include <fstream>
#include <string>
#include <map>
#include <set>

#include <pthread.h>

#define CFISH_USE_SHORT_NAMES
#define LUCY_USE_SHORT_NAMES
#include <Clownfish/String.h>

#include <Lucy/Analysis/RegexTokenizer.h>
#include <Lucy/Document/Doc.h>
#include <Lucy/Document/HitDoc.h>
#include <Lucy/Index/Indexer.h>
#include <Lucy/Plan/FullTextType.h>
#include <Lucy/Plan/Schema.h>
#include <Lucy/Plan/StringType.h>
#include <Lucy/Search/ANDQuery.h>
#include <Lucy/Search/Hits.h>
#include <Lucy/Search/IndexSearcher.h>
#include <Lucy/Search/QueryParser.h>
#include <Lucy/Search/TermQuery.h>

#include <GraphMol/AtomIterators.h>
#include <GraphMol/Fingerprints/IOCBFingerprint.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

/*
 * tools
 */

static const unsigned char b64str[65] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
static std::string fp2str (uint32_t fp)
{
	std::string r (6, '_');
	for (int i = 0; i < 6; ++i)
		r[i] = b64str[ (fp >> (6 * i)) & 0x3f];
	return r;
}

static Schema* create_schema (String*guidF, String*fpF)
{
	Schema*schema = Schema_new();

	{
		StringType *type = StringType_new();
		StringType_Set_Indexed (type, false);
		Schema_Spec_Field (schema, guidF, (FieldType*) type);
		DECREF (type);
	}

	{
		String* re = Str_newf ("[a-zA-Z0-9+/]+");
		RegexTokenizer *rt = RegexTokenizer_new (re);
		FullTextType *type = FullTextType_new ( (Analyzer*) rt);
		FullTextType_Set_Indexed (type, true);
		FullTextType_Set_Highlightable (type, false);
		FullTextType_Set_Stored (type, false);
		Schema_Spec_Field (schema, fpF, (FieldType*) type);
		DECREF (rt);
		DECREF (re);
		DECREF (type);
	}

	return schema;
}

/*
 * internal stuffs
 */

static std::map<std::string, int> fporder __attribute__ ((init_priority (150)));

void load_fporder()
{
	fporder.clear();
	std::ifstream sdf (DATADIR "/fporder.txt");
	std::string fp;
	int fpn = 0;
	while (std::getline (sdf, fp)) fporder[fp] = fpn++;
}

struct fpsearch_data {
	pthread_mutex_t buffer_mtx;
	std::list<Doc*> buffer;
	size_t bufsize;
	pthread_mutex_t index_mtx;

	String *folder, *guidF, *fpF, *boolop;
	Schema *schema;

	IndexSearcher *searcher;
	QueryParser *qparser;

	FPSEARCH_API (params_t) p;
};

static void commit_buffer (fpsearch_data&d, bool force)
{
	std::list<Doc*> batch;

	pthread_mutex_lock (&d.buffer_mtx);
	if (force || d.bufsize > d.p.indexerBatch) {
		batch.swap (d.buffer);
		d.bufsize = 0;
	}
	pthread_mutex_unlock (&d.buffer_mtx);

	if (batch.empty()) return;

	pthread_mutex_lock (&d.index_mtx);
	Indexer *indexer = Indexer_new (d.schema, (Obj*) (d.folder),
	                                nullptr, Indexer_CREATE);
	while (!batch.empty()) {
		Indexer_Add_Doc (indexer, batch.front(), 1.0);
		DECREF (batch.front());
		batch.pop_front();
	}
	Indexer_Commit (indexer);
	DECREF (indexer);
	pthread_mutex_unlock (&d.index_mtx);
}

static void index_optimize (fpsearch_data&d)
{
	commit_buffer (d, true);

	pthread_mutex_lock (&d.index_mtx);
	Indexer *indexer = Indexer_new (d.schema, (Obj*) (d.folder),
	                                nullptr, Indexer_CREATE);
	Indexer_Optimize (indexer);
	Indexer_Commit (indexer);
	DECREF (indexer);
	pthread_mutex_unlock (&d.index_mtx);
}

static void close_indexer (fpsearch_data&d)
{
	commit_buffer (d, true);
}

static void close_searcher (fpsearch_data&d)
{
	if (!d.searcher) return;
	pthread_mutex_lock (&d.index_mtx); //this may be called from add_mol
	if (d.searcher) {
		DECREF (d.searcher);
		d.searcher = nullptr;
	}
	pthread_mutex_unlock (&d.index_mtx);
}

static void open_indexer (fpsearch_data&d)
{
	close_searcher (d);
}

static void open_searcher (fpsearch_data&d)
{
	if (d.searcher) return; //already searching!
	close_indexer (d);
	d.searcher = IxSearcher_new ( (Obj*) (d.folder));
}

static inline RDKit::Bond::BondType bondtype_jg2rd (int multiplicity,
                                                    bool isAromatic)
{
	if (isAromatic) return RDKit::Bond::AROMATIC;
	switch (multiplicity) {
	case 1:
		return RDKit::Bond::SINGLE;
	case 2:
		return RDKit::Bond::DOUBLE;
	case 3:
		return RDKit::Bond::TRIPLE;
	case 4:
		return RDKit::Bond::QUADRUPLE;
	case 5:
		return RDKit::Bond::QUINTUPLE;
	case 6:
		return RDKit::Bond::HEXTUPLE;
	default:
		return RDKit::Bond::UNSPECIFIED;
	}
}

//caller responsible for deallocating the ROMol
static RDKit::ROMol* JGMol2RDMol (const Molecule*m)
{
	try {
		RDKit::RWMol rwm;
		std::map<int, int> atomMap;

		for (int i = 0; i < m->atomCount; ++i) {
			RDKit::Atom* a = new RDKit::Atom();
			int aNum = molecule_get_atom_number (m, i);
			a->setAtomicNum (aNum > 0 ? aNum : 0);
			a->setFormalCharge (molecule_get_formal_charge (m, i));
			a->setNumExplicitHs (molecule_get_hydrogen_count (m, i));
			a->setNoImplicit (true); //no implicit Hs
			atomMap[i] = rwm.addAtom (a, false, true);
		}

		for (int i = 0; i < m->atomCount; ++i) {
			int nbs = molecule_get_bond_list_size (m, i);
			for (int j = 0; j < nbs; ++j) {
				int otherAtom = molecule_get_bond_list (m, i) [j];
				if (otherAtom <= i) continue;
				int bondData = molecule_get_bond_data
				               (m, molecule_get_bond
				                (m, i, otherAtom));
				bool isAromatic = bondData & 0xc0;
				int multiplicity = (bondData & 0x7c) >> 3;
				rwm.addBond (atomMap[i], atomMap[otherAtom],
				             bondtype_jg2rd (multiplicity,
				                             isAromatic));
			}
		}

		RDKit::ROMol *rom = new RDKit::ROMol (rwm);
		return rom;
	} catch (...) {
		return nullptr;
	}
}

static void index_add (fpsearch_data&d, int guid, RDKit::ROMol*m)
{
	RDKit::SparseIntVect<uint32_t> *res
	    = RDKit::IOCBFingerprints::getFingerprint
	      (*m, d.p.graphSize, d.p.circSize, d.p.maxLogFeats);
	std::string fp;
	fp.reserve (res->size() * 7);
	for (auto&&i : res->getNonzeroElements())
		fp.append (fp2str (i.first) + " ");
	delete res;

	std::string guids = std::to_string (guid);

	pthread_mutex_lock (&d.index_mtx);
	Doc*doc = Doc_new (nullptr, 0);

	{
		String *value = Str_newf (guids.c_str());
		Doc_Store (doc, d.guidF, (Obj*) value);
		DECREF (value);
	}
	{
		String *value = Str_newf (fp.c_str());
		Doc_Store (doc, d.fpF, (Obj*) value);
		DECREF (value);
	}
	pthread_mutex_unlock (&d.index_mtx);

	pthread_mutex_lock (&d.buffer_mtx);
	d.buffer.push_back (doc);
	size_t bs = ++d.bufsize;
	pthread_mutex_unlock (&d.buffer_mtx);

	if (bs >= d.p.indexerBatch)
		commit_buffer (d, false);
}

static void index_remove (fpsearch_data&d, int guid)
{
	std::string guids = std::to_string (guid);
	pthread_mutex_lock (&d.index_mtx);

	Indexer *indexer = Indexer_new (d.schema, (Obj*) (d.folder),
	                                nullptr, Indexer_CREATE);
	String *value = Str_newf (guids.c_str());

	Indexer_Delete_By_Term (indexer, d.guidF, (Obj*) value);
	DECREF (value);

	Indexer_Commit (indexer);
	DECREF (indexer);

	pthread_mutex_unlock (&d.index_mtx);
}

static Hits* search_query (fpsearch_data&d, const Molecule*m, int max_results)
{
	auto* molh = JGMol2RDMol (m);
	if (!molh) return nullptr;
#if EXPLICITLY_REMOVE_HS
	auto *mol = RDKit::MolOps::removeHs (*molh, false, false, false);
	delete molh;
#else
	auto*mol = molh;
#endif

	RDKit::IOCBFingerprints::BitInfo bits;
	RDKit::SparseIntVect<uint32_t> *res =
	    RDKit::IOCBFingerprints::getFingerprint
	    (*mol, d.p.graphSize, d.p.circSize, d.p.maxLogFeats, true, &bits);
	int molAtoms = mol->getNumAtoms();
	delete mol;

	//convert and pre-sort the fingerprints
	std::map<int, std::pair<uint32_t, std::string> > fpi;
	int unknownId = -1;
	for (auto&&i : res->getNonzeroElements()) {
		std::string fp = fp2str (i.first);
		auto o = fporder.find (fp);
		if (o == fporder.end())
			/* this is a hack: if the fingerprint is not known to
			 * fporder (which it should be but keeping that
			 * database in shape isn't very easy), let's assume
			 * it's very good (and put it on the beginning of the
			 * queue... :] ) */
			fpi[unknownId--] = std::make_pair (i.first, fp);
		else
			fpi[o->second] = std::make_pair (i.first, fp);
	}
	delete res;

	//find a decent coverage
	std::set<std::string> fps;
	{
		std::vector<int> coverage;
		int uncovered = molAtoms, nfps = 0;
		coverage.resize (uncovered, 0);

		for (auto&i : fpi) {
			if (uncovered <= 0) break;
			if (nfps >= d.p.searchMaxFps) break;
			bool found = false;
			for (auto&a : bits[i.second.first]) {
				if (coverage[a] < d.p.searchAtomCoverage) {
					found = true;
					++coverage[a];
					if (coverage[a]
					    == d.p.searchAtomCoverage)
						--uncovered;
				}
			}
			if (found) {
				fps.insert (i.second.second);
				++nfps;
			}
		}
	}

	//convert the coverage to actual query
	std::string q;
	{
		bool first = true;
		q.reserve (fps.size() * 10);
		for (auto&i : fps) {
			if (!first) {
				q.append (" ");
			} else first = false;
			q.append ("\"" + i + "\"");
		}
	}

	//query -> lucy query
	String*query_str = Str_newf (q.c_str());
	Query *query = QParser_Parse (d.qparser, query_str);
	Hits*hits = IxSearcher_Hits (d.searcher, (Obj*) query, 0, max_results, nullptr);

	DECREF (query);
	DECREF (query_str);

	return hits;
}

/*
 * "interface"
 */

extern "C" {

	PG_MODULE_MAGIC;

	void __attribute__ ((constructor (160))) fpsearch_lucy_init (void)
	{
		load_fporder();
	}
}

void FPSEARCH_API (initialize) (void**ddp, const char*index_dir)
{
	lucy_bootstrap_parcel();

	*ddp = new fpsearch_data;
	fpsearch_data&d = * (fpsearch_data*) *ddp;

	d.buffer_mtx = PTHREAD_MUTEX_INITIALIZER;
	d.bufsize = 0;
	d.index_mtx = PTHREAD_MUTEX_INITIALIZER;

	d.folder = Str_newf (index_dir);
	d.guidF = Str_newf ("guid");
	d.fpF = Str_newf ("fp");
	d.boolop = Str_newf ("AND");
	d.schema = create_schema (d.guidF, d.fpF);
	d.searcher = nullptr;
	d.qparser = QParser_new (d.schema, nullptr, d.boolop, nullptr);

	d.p.graphSize = 7;
	d.p.circSize = 3;
	d.p.maxLogFeats = 5;
	d.p.searchAtomCoverage = 2;
	d.p.searchMaxFps = 32;
	d.p.indexerBatch = 10000;
}

void FPSEARCH_API (close) (void*dd)
{
	fpsearch_data&d = * (fpsearch_data*) dd;
	close_indexer (d);
	close_searcher (d);
	DECREF (d.folder);
	DECREF (d.guidF);
	DECREF (d.fpF);
	DECREF (d.boolop);
	DECREF (d.schema);
	DECREF (d.qparser);
	delete (fpsearch_data*) dd;
}

//retrieve the params
void FPSEARCH_API (params_get) (void*dd, FPSEARCH_API (params_t) *pp)
{
	fpsearch_data&d = * (fpsearch_data*) dd;
	*pp = d.p;
}

//set the params
void FPSEARCH_API (params_set) (void*dd, const FPSEARCH_API (params_t) *pp)
{
	fpsearch_data&d = * (fpsearch_data*) dd;
	d.p = *pp;
}

//returns a new search object
void* FPSEARCH_API (search) (void*dd, const Molecule*mol, int max_results)
{
	fpsearch_data&d = * (fpsearch_data*) dd;
	open_searcher (d);

	return (void*) search_query (d, mol, max_results);
}

//gets data from search object to array, max n, returns retrieved count
int FPSEARCH_API (search_fillbuf) (void*dd, void*ss, int*results, int n)
{
	fpsearch_data&d = * (fpsearch_data*) dd;
	Hits*hits = (Hits*) ss;
	int ret = 0;
	while (n > 0) {
		HitDoc*hit = Hits_Next (hits);
		if (!hit) break;
		String*guid = (String*) HitDoc_Extract (hit, d.guidF);
		* (results++) = Str_To_I64 (guid);
		++ret;
		--n;
		DECREF (guid);
		DECREF (hit);
	}
	return ret;
}

//remove the search result
void FPSEARCH_API (search_finish) (void*dd, void*ss)
{
	DECREF ( (Hits*) ss);
}

int FPSEARCH_API (add_mol) (void*dd, int guid, const Molecule*mol)
{
	fpsearch_data&d = * (fpsearch_data*) dd;
	open_indexer (d);
	auto*m = JGMol2RDMol (mol);
	if (!m) return 1;
	index_add (d, guid, m);
	delete m;
	return 0;
}

void FPSEARCH_API (remove_mol) (void*dd, int guid)
{
	fpsearch_data&d = * (fpsearch_data*) dd;
	open_indexer (d);
	index_remove (d, guid);
}

void FPSEARCH_API (flush) (void*dd)
{
	fpsearch_data&d = * (fpsearch_data*) dd;
	open_indexer (d);
	commit_buffer (d, true);
}

void FPSEARCH_API (optimize) (void*dd)
{
	fpsearch_data&d = * (fpsearch_data*) dd;
	open_indexer (d);
	index_optimize (d);
}
