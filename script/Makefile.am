bin_SCRIPTS = \
	sachem-updatedb

EXTRA_DIST = \
	sachem-updatedb.in

CLEANFILES = $(bin_SCRIPTS)


sachem-updatedb: sachem-updatedb.in Makefile
		sed "s|@CLASSPATH@|$(datadir)/$(PACKAGE)/jar|" < $(srcdir)/sachem-updatedb.in > sachem-updatedb
		chmod +x sachem-updatedb
