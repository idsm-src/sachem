pglibdir = $(POSTGRESQL_LIBDIR)
pglib_LTLIBRARIES = libsachem.la

libsachem_la_CFLAGS = -std=gnu99 -O3 -march=native -flto -minline-stringops-dynamically -Wall -g
libsachem_la_CXXFLAGS = -std=c++11 -O3 -march=native -flto -minline-stringops-dynamically -Wall -g
libsachem_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS) $(JDK8_CPPFLAGS) -DBINDIR=\"$(bindir)\" -DJARDIR=\"$(datadir)/$(PACKAGE)/jar\"
libsachem_la_LDFLAGS = $(JDK8_LDFLAGS) -Wl,-rpath,$(POSTGRESQL_LIBDIR)

libsachem_la_SOURCES = \
		info.c \
        java.c \
        native.c \
        sachem.c \
        search.c \
        sync.c

EXTRA_DIST = \
        enum.h \
        isomorphism.h \
        java.h \
        memory.h \
        molecule.h \
        native.h \
        sachem.h
