bin_SCRIPTS = inchi/inchi-1

jar_DATA = lib/cdk-2.1.1.jar sachem.jar lib/commons-net-3.6.jar lib/jdbc-postgresql.jar lib/lucene-core-7.2.1.jar lib/lucene-queries-7.2.1.jar
jardir = $(datadir)/$(PACKAGE)/jar

sachem.jar: force-build

force-build:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) build-jar

clean-local:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) clean

EXTRA_DIST = \
	build.xml \
	inchi/inchi-1 \
    lib/cdk-2.1.1.jar \
    lib/commons-net-3.6.jar \
    lib/jdbc-postgresql.jar \
    lib/lucene-core-7.2.1.jar \
    lib/lucene-queries-7.2.1.jar \
    src/cz/iocb/sachem/convert/ConvertMolecule.java \
    src/cz/iocb/sachem/fingerprint/bitpos/BitPosApi.java \
    src/cz/iocb/sachem/fingerprint/bitpos/BitPositions.java \
    src/cz/iocb/sachem/fingerprint/bitpos/ExtendedBitPositions.java \
    src/cz/iocb/sachem/fingerprint/bitpos/Neighbour.java \
    src/cz/iocb/sachem/fingerprint/ExtendedFingerprinter.java \
    src/cz/iocb/sachem/fingerprint/OrchemExtendedFingerprinter.java \
    src/cz/iocb/sachem/fingerprint/OrchemFingerprinter.java \
    src/cz/iocb/sachem/isomorphism/IsomorphismSort.java \
    src/cz/iocb/sachem/load/ConfigurationProperties.java \
    src/cz/iocb/sachem/load/CompoundLoader.java \
    src/cz/iocb/sachem/load/CompoundUpdater.java \
    src/cz/iocb/sachem/load/ConnectionPool.java \
    src/cz/iocb/sachem/load/GeneralCompoundLoader.java \
    src/cz/iocb/sachem/load/DrugbankCompoundUpdater.java \
    src/cz/iocb/sachem/load/PubChemCompoundUpdater.java \
    src/cz/iocb/sachem/lucene/BitSetCollector.java \
    src/cz/iocb/sachem/lucene/FingerprintAnalyzer.java \
    src/cz/iocb/sachem/lucene/FingerprintReader.java \
    src/cz/iocb/sachem/lucene/FingerprintTokenizer.java \
    src/cz/iocb/sachem/lucene/Lucene.java \
    src/cz/iocb/sachem/lucene/ScoreHit.java \
    src/cz/iocb/sachem/lucene/SimilarDocCollector.java \
    src/cz/iocb/sachem/lucene/TanimotoScoreProvider.java \
    src/cz/iocb/sachem/lucene/TanimotoScoreQuery.java \
    src/cz/iocb/sachem/search/LucyLoader.java \
    src/cz/iocb/sachem/search/EcdkLoader.java \
    src/cz/iocb/sachem/search/EcdkSimilaritySearch.java \
    src/cz/iocb/sachem/search/EcdkSubstructureSearch.java \
    src/cz/iocb/sachem/search/OrchemLoader.java \
    src/cz/iocb/sachem/search/OrchemSimilaritySearch.java \
    src/cz/iocb/sachem/search/OrchemSubstructureSearch.java \
    src/cz/iocb/sachem/search/SachemMoleculeBuilder.java \
    src/cz/iocb/sachem/search/SimilaritySearch.java \
    src/cz/iocb/sachem/search/SubstructureSearch.java \
    src/cz/iocb/sachem/search/QueryFormat.java \
    src/cz/iocb/sachem/shared/MoleculeCounts.java \
    src/cz/iocb/sachem/shared/MoleculeCreator.java \
    src/cz/iocb/sachem/tautomers/CombinationCountException.java \
    src/cz/iocb/sachem/tautomers/InChI.java \
    src/cz/iocb/sachem/tautomers/InchiTautomerGenerator.java
