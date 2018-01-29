pglibdir = $(POSTGRESQL_LIBDIR)
pglib_LTLIBRARIES = libfingerprints.la

libfingerprints_la_CPPFLAGS = -I$(srcdir)/.. $(POSTGRESQL_CPPFLAGS) -DDATADIR=\"$(datadir)/$(PACKAGE)\"
libfingerprints_la_CXXFLAGS = -std=c++11 -Ofast -march=native

libfingerprints_la_SOURCES = \
	fingerprint.cpp \
	IOCBFingerprint.cpp \
	AtomFingerprint.cpp \
	CRNGFingerprint.cpp \
	SGFingerprint.cpp

EXTRA_DIST = \
	fingerprint.h \
	IOCBFingerprint.hpp \
	AtomFingerprint.hpp \
	CRNGFingerprint.hpp \
	SGFingerprint.hpp \
	SubstructureMatch.hpp
