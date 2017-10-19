lib_LTLIBRARIES = libfpsearch.la

libfpsearch_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(CLOWNFISH_CPPFLAGS) $(LUCY_CPPFLAGS) $(RDKIT_CPPFLAGS)
libfpsearch_la_CXXFLAGS = -std=c++11
libfpsearch_la_LIBADD = -llucy -lclownfish
libfpsearch_la_LDFLAGS = $(CLOWNFISH_LDFLAGS) $(LUCY_LDFLAGS) \
		../../build/rdkit/Code/DataStructs/libRDKitDataStructs_static.a \
		../../build/rdkit/Code/GraphMol/Fingerprints/libRDKitFingerprints_static.a \
		../../build/rdkit/Code/GraphMol/Subgraphs/libRDKitSubgraphs_static.a \
		../../build/rdkit/Code/GraphMol/Substruct/libRDKitSubstructMatch_static.a \
		../../build/rdkit/Code/GraphMol/libRDKitGraphMol_static.a \
		../../build/rdkit/Code/GraphMol/SmilesParse/libRDKitSmilesParse_static.a \
		../../build/rdkit/Code/Geometry/libRDKitRDGeometryLib_static.a \
		../../build/rdkit/Code/RDGeneral/libRDKitRDGeneral_static.a

libfpsearch_la_SOURCES = fpsearch-lucy.cpp
EXTRA_DIST = fpsearch-lucy.h
