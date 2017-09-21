bin_SCRIPTS = inchi/inchi-1

jar_DATA = lib/cdk-2.0.jar orchem.jar
jardir = $(prefix)/$(PACKAGE)/jar

orchem.jar: force-build

force-build:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) build-jar

clean-local:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) clean

EXTRA_DIST = \
	build.xml \
	lib \
	src
