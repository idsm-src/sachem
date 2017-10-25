lib_LTLIBRARIES = libfpsearch.la

libfpsearch_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(CLOWNFISH_CPPFLAGS) $(LUCY_CPPFLAGS) $(RDKIT_CPPFLAGS)
libfpsearch_la_CXXFLAGS = -std=c++11
libfpsearch_la_LIBADD = -llucy -lclownfish -lboost_thread
libfpsearch_la_LDFLAGS = $(CLOWNFISH_LDFLAGS) $(LUCY_LDFLAGS) \
		../../build/rdkit/Code/DataStructs/CMakeFiles/DataStructs.dir/*.o \
		../../build/rdkit/Code/GraphMol/Fingerprints/CMakeFiles/Fingerprints.dir/*.o \
		../../build/rdkit/Code/GraphMol/Subgraphs/CMakeFiles/Subgraphs.dir/*.o \
		../../build/rdkit/Code/GraphMol/Substruct/CMakeFiles/SubstructMatch.dir/*.o \
		../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/*.o \
		../../build/rdkit/Code/GraphMol/SmilesParse/CMakeFiles/SmilesParse.dir/*.o \
		../../build/rdkit/Code/Geometry/CMakeFiles/RDGeometryLib.dir/*.o \
		../../build/rdkit/Code/RDGeneral/CMakeFiles/RDGeneral.dir/*.o

libfpsearch_la_SOURCES = fpsearch-lucy.cpp
EXTRA_DIST = fpsearch-lucy.h
