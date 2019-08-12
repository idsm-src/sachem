pglibdir = $(POSTGRESQL_LIBDIR)
pglib_LTLIBRARIES = libsachem.la

libsachem_la_CFLAGS = -std=gnu99 -O3 -march=native -flto -minline-stringops-dynamically -Wall -g
libsachem_la_CXXFLAGS = -std=c++11 -O3 -march=native -flto -minline-stringops-dynamically -Wall -g
libsachem_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(JDK8_CPPFLAGS) $(CLOWNFISH_CPPFLAGS) $(LUCY_CPPFLAGS) -DBINDIR=\"$(bindir)\" -DJARDIR=\"$(datadir)/$(PACKAGE)/jar\"
libsachem_la_LDFLAGS = $(JDK8_LDFLAGS) $(CLOWNFISH_LDFLAGS) $(LUCY_LDFLAGS) -Wl,-rpath,$(POSTGRESQL_LIBDIR)
libsachem_la_LIBADD = -llucy -lclownfish

libsachem_la_SOURCES = \
		isomorphism.c \
		fporder.c \
        molindex.c \
        sachem.c \
        stats.cpp \
        fingerprints/fingerprint.cpp \
	    fingerprints/IOCBFingerprint.cpp \
	    fingerprints/AtomFingerprint.cpp \
	    fingerprints/CRNGFingerprint.cpp \
	    fingerprints/RCFingerprint.cpp \
	    fingerprints/SGFingerprint.cpp \
        java/java.c \
        java/parse.c \
        ecdk/ecdk.c \
        ecdk/simsearch.c \
        ecdk/subsearch.c \
        ecdk/sync.c \
        lucene/indexer.c \
        lucene/lucene.c \
        lucene/search.c \
        lucene/simsearch.c \
        lucene/subsearch.c \
        lucene/sync.c \
        lucy/lucy.c \
        lucy/subsearch.c \
        lucy/sync.c \
        orchem/orchem.c \
        orchem/simsearch.c \
        orchem/subsearch.c \
        orchem/sync.c

EXTRA_DIST = \
        bitset.h \
        heap.h \
        fporder.h \
        isomorphism.h \
        measurement.h \
        molecule.h \
        molindex.h \
        sachem.h \
        stats.h \
        subsearch.h \
        fingerprints/fingerprint.h \
	    fingerprints/IOCBFingerprint.hpp \
	    fingerprints/AtomFingerprint.hpp \
	    fingerprints/CRNGFingerprint.hpp \
	    fingerprints/RCFingerprint.hpp \
	    fingerprints/SGFingerprint.hpp \
	    fingerprints/SubstructureMatch.hpp \
        java/java.h \
        java/parse.h \
        ecdk/common.h \
        ecdk/ecdk.h \
        lucene/common.h \
        lucene/indexer.h \
        lucene/lucene.h \
        lucene/search.h \
        lucy/common.h \
        lucy/lucy.h \
        orchem/common.h \
        orchem/orchem.h
