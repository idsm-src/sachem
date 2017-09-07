lib_LTLIBRARIES = libpgchem.la

libpgchem_la_CPPFLAGS = -std=gnu99 -Ofast -march=native $(POSTGRESQL_CPPFLAGS) $(JDK8_CPPFLAGS) -DJARDIR=\"$(prefix)/$(PACKAGE)/jar\"
libpgchem_la_LDFLAGS = $(JDK8_LDFLAGS)

libpgchem_la_SOURCES = \
		pgchem.c \
		java.c \
		orchem.c
