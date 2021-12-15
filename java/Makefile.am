jar_DATA = \
	sachem.jar \
    lib/beam-core-1.3.3.jar \
    lib/cdk-atomtype-2.3.jar \
    lib/cdk-core-2.3.jar \
    lib/cdk-ctab-2.3.jar \
    lib/cdk-interfaces-2.3.jar \
    lib/cdk-ioformats-2.3.jar \
    lib/cdk-isomorphism-2.3.jar \
    lib/cdk-silent-2.3.jar \
    lib/cdk-smiles-2.3.jar \
    lib/cdk-standard-2.3.jar \
    lib/cdk-valencycheck-2.3.jar \
    lib/commons-net-3.6.jar \
    lib/guava-29.0-jre.jar \
    lib/log4j-1.2-api-2.16.0.jar \
    lib/log4j-api-2.16.0.jar \
    lib/log4j-core-2.16.0.jar \
    lib/lucene-core-8.7.0.jar \
    lib/postgresql-42.2.14.jar \
    lib/vecmath-1.5.2.jar

jardir = $(datadir)/$(PACKAGE)/jar

sachem.jar: force-build

force-build:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) build-jar

clean-local:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) clean

EXTRA_DIST = \
	build.xml \
	patterns.bin \
    lib/beam-core-1.3.3.jar \
    lib/cdk-atomtype-2.3.jar \
    lib/cdk-core-2.3.jar \
    lib/cdk-ctab-2.3.jar \
    lib/cdk-interfaces-2.3.jar \
    lib/cdk-ioformats-2.3.jar \
    lib/cdk-isomorphism-2.3.jar \
    lib/cdk-silent-2.3.jar \
    lib/cdk-smiles-2.3.jar \
    lib/cdk-standard-2.3.jar \
    lib/cdk-valencycheck-2.3.jar \
    lib/commons-net-3.6.jar \
    lib/guava-29.0-jre.jar \
    lib/log4j-1.2-api-2.16.0.jar \
    lib/log4j-api-2.16.0.jar \
    lib/log4j-core-2.16.0.jar \
    lib/lucene-core-8.7.0.jar \
    lib/postgresql-42.2.14.jar \
    lib/vecmath-1.5.2.jar \
    src/cz/iocb/sachem/fingerprint/AtomFingerprint.java \
    src/cz/iocb/sachem/fingerprint/CRNGFingerprint.java \
    src/cz/iocb/sachem/fingerprint/Fingerprint.java \
    src/cz/iocb/sachem/fingerprint/IOCBFingerprint.java \
    src/cz/iocb/sachem/fingerprint/RCFingerprint.java \
    src/cz/iocb/sachem/fingerprint/SGFingerprint.java \
    src/cz/iocb/sachem/load/CompoundLoader.java \
    src/cz/iocb/sachem/load/CompoundUpdater.java \
    src/cz/iocb/sachem/load/ConfigurationProperties.java \
    src/cz/iocb/sachem/load/DrugbankCompoundUpdater.java \
    src/cz/iocb/sachem/load/GeneralCompoundLoader.java \
    src/cz/iocb/sachem/load/PubChemCompoundUpdater.java \
    src/cz/iocb/sachem/lucene/BalancedMergePolicy.java \
    src/cz/iocb/sachem/lucene/FingerprintBitMapping.java \
    src/cz/iocb/sachem/lucene/FingerprintTokenStream.java \
    src/cz/iocb/sachem/lucene/Indexer.java \
    src/cz/iocb/sachem/lucene/IndexInfo.java \
    src/cz/iocb/sachem/lucene/ResultCollectorManager.java \
    src/cz/iocb/sachem/lucene/Searcher.java \
    src/cz/iocb/sachem/lucene/SearchResult.java \
    src/cz/iocb/sachem/lucene/Settings.java \
    src/cz/iocb/sachem/lucene/SimilarStructureQuery.java \
    src/cz/iocb/sachem/lucene/SortedResultCollectorManager.java \
    src/cz/iocb/sachem/lucene/SubstructureQuery.java \
    src/cz/iocb/sachem/lucene/TopResultCollectorManager.java \
    src/cz/iocb/sachem/molecule/AromaticityMode.java \
    src/cz/iocb/sachem/molecule/BinaryMoleculeBuilder.java \
    src/cz/iocb/sachem/molecule/BinaryMolecule.java \
    src/cz/iocb/sachem/molecule/BinaryMoleculeSort.java \
    src/cz/iocb/sachem/molecule/ChargeMode.java \
    src/cz/iocb/sachem/molecule/InChITautomerGenerator.java \
    src/cz/iocb/sachem/molecule/InChITools.java \
    src/cz/iocb/sachem/molecule/Isomorphism.java \
    src/cz/iocb/sachem/molecule/IsotopeMode.java \
    src/cz/iocb/sachem/molecule/MoleculeCreator.java \
    src/cz/iocb/sachem/molecule/Molecule.java \
    src/cz/iocb/sachem/molecule/NativeIsomorphism.java \
    src/cz/iocb/sachem/molecule/QueryFormat.java \
    src/cz/iocb/sachem/molecule/RadicalMode.java \
    src/cz/iocb/sachem/molecule/SearchMode.java \
    src/cz/iocb/sachem/molecule/StereoMode.java \
    src/cz/iocb/sachem/molecule/TautomerMode.java
