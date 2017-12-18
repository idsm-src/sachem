bin_SCRIPTS = \
	sachem-updatedb \
	sachem-update-drugbank

EXTRA_DIST = \
	sachem-updatedb.in \
	sachem-update-drugbank.in

CLEANFILES = $(bin_SCRIPTS)


sachem-%: sachem-%.in Makefile
		sed "s|@CLASSPATH@|$(datadir)/$(PACKAGE)/jar|" < $< > $@
		chmod +x $@
