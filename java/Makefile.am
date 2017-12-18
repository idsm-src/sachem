bin_SCRIPTS = inchi/inchi-1

jar_DATA = lib/cdk-2.0.jar sachem.jar lib/commons-net-3.6.jar lib/jdbc-postgresql.jar
jardir = $(datadir)/$(PACKAGE)/jar

sachem.jar: force-build

force-build:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) build-jar

clean-local:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) clean

EXTRA_DIST = \
	build.xml \
	inchi/inchi-1 \
    lib/cdk-2.0.jar \
    lib/commons-net-3.6.jar \
    lib/jdbc-postgresql.jar \
    src/cz/iocb/orchem/convert/ConvertMolecule.java \
    src/cz/iocb/orchem/fingerprint/bitpos/BitPosApi.java \
    src/cz/iocb/orchem/fingerprint/bitpos/BitPositions.java \
    src/cz/iocb/orchem/fingerprint/bitpos/ExtendedBitPositions.java \
    src/cz/iocb/orchem/fingerprint/bitpos/Neighbour.java \
    src/cz/iocb/orchem/fingerprint/OrchemExtendedFingerprinter.java \
    src/cz/iocb/orchem/fingerprint/OrchemFingerprinter.java \
    src/cz/iocb/orchem/isomorphism/IsomorphismSort.java \
    src/cz/iocb/orchem/load/ConfigurationProperties.java \
    src/cz/iocb/orchem/load/CompoundLoader.java \
    src/cz/iocb/orchem/load/CompoundUpdater.java \
    src/cz/iocb/orchem/load/ConnectionPool.java \
    src/cz/iocb/orchem/load/ChebiCompoundLoader.java \
    src/cz/iocb/orchem/load/ChemblCompoundLoader.java \
    src/cz/iocb/orchem/load/DrugbankCompoundUpdater.java \
    src/cz/iocb/orchem/load/PubChemCompoundLoader.java \
    src/cz/iocb/orchem/search/OrchemLoader.java \
    src/cz/iocb/orchem/search/OrchemMoleculeBuilder.java \
    src/cz/iocb/orchem/search/OrchemSimilaritySearch.java \
    src/cz/iocb/orchem/search/OrchemSubstructureSearch.java \
    src/cz/iocb/orchem/search/SubstructureSearch.java \
    src/cz/iocb/orchem/shared/MoleculeCounts.java \
    src/cz/iocb/orchem/shared/MoleculeCreator.java \
    src/cz/iocb/orchem/tautomers/CombinationCountException.java \
    src/cz/iocb/orchem/tautomers/InChI.java \
    src/cz/iocb/orchem/tautomers/InchiTautomerGenerator.java
