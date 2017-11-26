lib_LTLIBRARIES = libsachem.la

libsachem_la_CPPFLAGS = -std=gnu99 -Ofast -march=native $(POSTGRESQL_CPPFLAGS) $(JDK8_CPPFLAGS) -DBINDIR=\"$(bindir)\" -DJARDIR=\"$(datadir)/$(PACKAGE)/jar\"
libsachem_la_LDFLAGS = $(JDK8_LDFLAGS)

libsachem_la_SOURCES = \
		isomorphism.c \
		java.c \
		sachem.c \
		simsearch.c \
		subsearch.c \
		subsearch-gin.c

EXTRA_DIST = \
		bitset.h \
		heap.h \
        isomorphism.h \
        java.h \
        molecule.h \
        sachem.h \
        simsearch.h \
        subsearch.h \
        subsearch-gin.h
