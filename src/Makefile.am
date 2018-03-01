SUBDIRS=fingerprints

pglibdir = $(POSTGRESQL_LIBDIR)
pglib_LTLIBRARIES = libsachem.la

libsachem_la_CFLAGS = -std=gnu99 -Ofast -march=native
libsachem_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(JDK8_CPPFLAGS) $(CLOWNFISH_CPPFLAGS) $(LUCY_CPPFLAGS) -DBINDIR=\"$(bindir)\" -DDATADIR=\"$(datadir)/$(PACKAGE)\" -DJARDIR=\"$(datadir)/$(PACKAGE)/jar\"
libsachem_la_LDFLAGS = $(JDK8_LDFLAGS) $(CLOWNFISH_LDFLAGS) $(LUCY_LDFLAGS)
libsachem_la_LIBADD = fingerprints/libfingerprints.la -llucy -lclownfish

libsachem_la_SOURCES = \
		isomorphism.c \
		fporder.c \
        java.c \
        lucene.c \
        lucy.c \
        sachem.c \
        simsearch.c \
        stats.cpp \
        subsearch.c \
        subsearch-lucene.c \
        subsearch-lucy.c

EXTRA_DIST = \
        bitset.h \
        heap.h \
        isomorphism.h \
        java.h \
        lucene.h \
        lucy.h \
        molecule.h \
        sachem.h \
        stats.h \
        subsearch.h
