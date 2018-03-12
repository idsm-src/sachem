pglibdir = $(POSTGRESQL_LIBDIR)
pglib_LTLIBRARIES = libsachem.la

libsachem_la_CFLAGS = -std=gnu99 -O3 -march=native -flto -minline-stringops-dynamically -Wall -g
libsachem_la_CXXFLAGS = -std=c++11 -O3 -march=native -flto -minline-stringops-dynamically -Wall -g
libsachem_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(JDK8_CPPFLAGS) $(CLOWNFISH_CPPFLAGS) $(LUCY_CPPFLAGS) -DBINDIR=\"$(bindir)\" -DJARDIR=\"$(datadir)/$(PACKAGE)/jar\"
libsachem_la_LDFLAGS = $(JDK8_LDFLAGS) $(CLOWNFISH_LDFLAGS) $(LUCY_LDFLAGS)
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
	    fingerprints/SGFingerprint.cpp \
        java/java.c \
        java/parse.c \
        lucene/lucene.c \
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
	    fingerprints/SGFingerprint.hpp \
	    fingerprints/SubstructureMatch.hpp \
        java/java.h \
        java/parse.h \
        lucene/common.h \
        lucene/lucene.h \
        lucy/common.h \
        lucy/lucy.h \
        orchem/common.h \
        orchem/orchem.h
