pglibdir = $(POSTGRESQL_LIBDIR)
pglib_LTLIBRARIES = libsachem.la

libsachem_la_CFLAGS = -std=gnu99 -O3 -flto -Wall -Wno-parentheses -g
libsachem_la_CXXFLAGS = -std=c++11 -O3 -flto -Wall -g
libsachem_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(JDK8_CPPFLAGS) -I$(srcdir)/inchi -DTARGET_API_LIB -DBINDIR=\"$(bindir)\" -DJARDIR=\"$(datadir)/$(PACKAGE)/jar\"
libsachem_la_LDFLAGS = $(JDK8_LDFLAGS) -Wl,-rpath,$(POSTGRESQL_LIBDIR)

libsachem_la_SOURCES = \
        inchi/ichi_bns.c \
        inchi/ichicano.c \
        inchi/ichicans.c \
        inchi/ichican2.c \
        inchi/ichiisot.c \
        inchi/ichimake.c \
        inchi/ichimap1.c \
        inchi/ichimap2.c \
        inchi/ichimap4.c \
        inchi/ichinorm.c \
        inchi/ichiqueu.c \
        inchi/ichiring.c \
        inchi/ichisort.c \
        inchi/ichister.c \
        inchi/ichitaut.c \
        inchi/readinch.c \
        inchi/strutil.c \
        inchi/util.c \
        inchi.c \
        info.c \
        isomorphism.c \
        java.c \
        sachem.c \
        search.c \
        sync.c

EXTRA_DIST = \
        inchi/extr_ct.h \
        inchi/ichi_bns.h \
        inchi/ichicano.h \
        inchi/ichicant.h \
        inchi/ichicomn.h \
        inchi/ichicomp.h \
        inchi/ichidrp.h \
        inchi/ichierr.h \
        inchi/ichi.h \
        inchi/ichinorm.h \
        inchi/ichiring.h \
        inchi/ichisize.h \
        inchi/ichister.h \
        inchi/ichitaut.h \
        inchi/ichitime.h \
        inchi/incomdef.h \
        inchi/inchi_api.h \
        inchi/inpdef.h \
        inchi/mode.h \
        inchi/readinch.h \
        inchi/strutil.h \
        inchi/util.h \
        enum.h \
        inchi.h \
        isomorphism.h \
        java.h \
        memory.h \
        molecule.h \
        sachem.h
