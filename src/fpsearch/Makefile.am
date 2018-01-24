pglibdir = $(POSTGRESQL_LIBDIR)
pglib_LTLIBRARIES = libfpsearch.la

libfpsearch_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(CLOWNFISH_CPPFLAGS) $(LUCY_CPPFLAGS) $(RDKIT_CPPFLAGS) -DDATADIR=\"$(datadir)/$(PACKAGE)\"

libfpsearch_la_CXXFLAGS = -std=c++11
libfpsearch_la_LIBADD = -llucy -lclownfish -lboost_thread
libfpsearch_la_LDFLAGS = $(CLOWNFISH_LDFLAGS) $(LUCY_LDFLAGS) \
        ../../build/rdkit/Code/Geometry/CMakeFiles/RDGeometryLib.dir/point.cpp.o \
        ../../build/rdkit/Code/Geometry/CMakeFiles/RDGeometryLib.dir/Transform3D.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/AddHs.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/Aromaticity.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/Atom.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/atomic_data.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/AtomIterators.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/Bond.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/BondIterators.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/Conformer.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/ConjugHybrid.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/FindRings.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/Chirality.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/Kekulize.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/Matrices.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/MolOps.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/MolPickler.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/PeriodicTable.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/QueryAtom.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/QueryBond.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/QueryOps.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/Resonance.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/RingInfo.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/ROMol.cpp.o \
        ../../build/rdkit/Code/GraphMol/CMakeFiles/GraphMol.dir/RWMol.cpp.o \
        ../../build/rdkit/Code/GraphMol/SmilesParse/CMakeFiles/SmilesParse.dir/CXSmilesOps.cpp.o \
        ../../build/rdkit/Code/GraphMol/SmilesParse/CMakeFiles/SmilesParse.dir/lex.yysmarts.cpp.o \
        ../../build/rdkit/Code/GraphMol/SmilesParse/CMakeFiles/SmilesParse.dir/lex.yysmiles.cpp.o \
        ../../build/rdkit/Code/GraphMol/SmilesParse/CMakeFiles/SmilesParse.dir/smarts.tab.cpp.o \
        ../../build/rdkit/Code/GraphMol/SmilesParse/CMakeFiles/SmilesParse.dir/SmilesParse.cpp.o \
        ../../build/rdkit/Code/GraphMol/SmilesParse/CMakeFiles/SmilesParse.dir/SmilesParseOps.cpp.o \
        ../../build/rdkit/Code/GraphMol/SmilesParse/CMakeFiles/SmilesParse.dir/smiles.tab.cpp.o \
        ../../build/rdkit/Code/GraphMol/Subgraphs/CMakeFiles/Subgraphs.dir/Subgraphs.cpp.o \
        ../../build/rdkit/Code/GraphMol/Subgraphs/CMakeFiles/Subgraphs.dir/SubgraphUtils.cpp.o \
        ../../build/rdkit/Code/GraphMol/Substruct/CMakeFiles/SubstructMatch.dir/SubstructMatch.cpp.o \
        ../../build/rdkit/Code/GraphMol/Substruct/CMakeFiles/SubstructMatch.dir/SubstructUtils.cpp.o \
        ../../build/rdkit/Code/RDGeneral/CMakeFiles/RDGeneral.dir/Dict.cpp.o \
        ../../build/rdkit/Code/RDGeneral/CMakeFiles/RDGeneral.dir/Invariant.cpp.o \
        ../../build/rdkit/Code/RDGeneral/CMakeFiles/RDGeneral.dir/LocaleSwitcher.cpp.o \
        ../../build/rdkit/Code/RDGeneral/CMakeFiles/RDGeneral.dir/RDLog.cpp.o \
        ../../build/rdkit/Code/RDGeneral/CMakeFiles/RDGeneral.dir/types.cpp.o \
        ../../build/rdkit/Code/RDGeneral/CMakeFiles/RDGeneral.dir/utils.cpp.o \
        ../../build/rdkit/Code/RDGeneral/CMakeFiles/RDGeneral.dir/versions.cpp.o

libfpsearch_la_SOURCES = \
	fpsearch-lucy.cpp \
	fingerprints/IOCBFingerprint.cpp \
	fingerprints/AtomFingerprint.cpp \
	fingerprints/CRNGFingerprint.cpp \
	fingerprints/SGFingerprint.cpp

EXTRA_DIST = \
	fpsearch-lucy.h \
	fingerprints/IOCBFingerprint.h \
	fingerprints/AtomFingerprint.h \
	fingerprints/CRNGFingerprint.h \
	fingerprints/SGFingerprint.h
