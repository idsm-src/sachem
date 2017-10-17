lib_LTLIBRARIES = libfpsearch.la

libfpsearch_la_CPPFLAGS = $(LUCY_CPPFLAGS) $(RDKIT_CPPFLAGS)
libfpsearch_la_CXXFLAGS = -std=c++11
libfpsearch_la_LDFLAGS = $(LUCY_LDFLAGS) $(RDKIT_LDFLAGS)
libfpsearch_la_LDADD = $(LUCY_LDADD) $(RDKIT_LDADD)

libfpsearch_la_SOURCES = fpsearch-lucy.cpp
EXTRA_DIST = fpsearch-lucy.h


