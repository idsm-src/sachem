pglibdir = $(POSTGRESQL_LIBDIR)
pglib_LTLIBRARIES = libsachem.la

libsachem_la_CFLAGS = -std=gnu99 -O3 -flto -Wall -g
libsachem_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(JDK8_CPPFLAGS) -I$(srcdir)/../lib/inchi/INCHI-1-SRC/INCHI_BASE/src -DTARGET_API_LIB -DCOMPILE_ANSI_ONLY -DBINDIR=\"$(bindir)\" -DJARDIR=\"$(datadir)/$(PACKAGE)/jar\"
libsachem_la_LDFLAGS = $(JDK8_LDFLAGS) -Wl,-rpath,$(POSTGRESQL_LIBDIR)
libsachem_la_LIBADD = ../lib/libinchi.la

libsachem_la_SOURCES = \
        inchi.c \
        info.c \
        isomorphism.c \
        java.c \
        sachem.c \
        search.c \
        sync.c

EXTRA_DIST = \
        enum.h \
        inchi.h \
        isomorphism.h \
        java.h \
        memory.h \
        molecule.h \
        sachem.h
