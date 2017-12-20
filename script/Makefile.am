bin_SCRIPTS = \
	sachem-update-data \
	sachem-update-drugbank \
	sachem-update-pubchem

EXTRA_DIST = \
	sachem-update-data.in \
	sachem-update-drugbank.in \
	sachem-update-pubchem.in

CLEANFILES = $(bin_SCRIPTS)


sachem-%: sachem-%.in Makefile
		sed "s|@CLASSPATH@|$(datadir)/$(PACKAGE)/jar|" < $< > $@
		chmod +x $@
