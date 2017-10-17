SUBDIRS=fpsearch

lib_LTLIBRARIES = libpgchem.la

libpgchem_la_CFLAGS = -std=gnu99 -Ofast -march=native
libpgchem_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(JDK8_CPPFLAGS) -DBINDIR=\"$(bindir)\" -DJARDIR=\"$(prefix)/$(PACKAGE)/jar\"
libpgchem_la_LDFLAGS = $(JDK8_LDFLAGS)
libpgchem_la_LIBADD = fpsearch/libfpsearch.la

libpgchem_la_SOURCES = \
        java.c \
        pgchem.c \
        simsearch.c \
        subsearch.c \
        subsearch-gin.c

EXTRA_DIST = \
        bitset.h \
        heap.h \
        isomorphism.h \
        java.h \
        molecule.h \
        pgchem.h \
        simsearch.h \
        subsearch.h \
        subsearch-gin.h
