jar_DATA = \
	sachem.jar \
    lib/beam-core-1.3.3.jar \
    lib/cdk-atomtype-2.3.jar \
    lib/cdk-core-2.3.jar \
    lib/cdk-ctab-2.3.jar \
    lib/cdk-inchi-2.3.jar \
    lib/cdk-interfaces-2.3.jar \
    lib/cdk-ioformats-2.3.jar \
    lib/cdk-isomorphism-2.3.jar \
    lib/cdk-silent-2.3.jar \
    lib/cdk-smiles-2.3.jar \
    lib/cdk-standard-2.3.jar \
    lib/cdk-valencycheck-2.3.jar \
    lib/commons-net-3.6.jar \
    lib/guava-28.1-jre.jar \
    lib/jdbc-postgresql.jar \
    lib/jnati-core-0.4.jar \
    lib/jnati-deploy-0.4.jar \
    lib/jni-inchi-0.8.jar \
    lib/log4j-1.2.17.jar \
    lib/lucene-core-8.2.0.jar \
    lib/vecmath-1.5.2.jar

jardir = $(datadir)/$(PACKAGE)/jar

sachem.jar: force-build

force-build:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) build-jar

clean-local:
		$(ANT) -f $(srcdir)/build.xml -Dbasedir=$(builddir) -Dbuild.dir=$(builddir) -Dsrc.dir=$(srcdir) clean

EXTRA_DIST = \
	build.xml \
    lib/beam-core-1.3.3.jar \
    lib/cdk-atomtype-2.3.jar \
    lib/cdk-core-2.3.jar \
    lib/cdk-ctab-2.3.jar \
    lib/cdk-inchi-2.3.jar \
    lib/cdk-interfaces-2.3.jar \
    lib/cdk-ioformats-2.3.jar \
    lib/cdk-isomorphism-2.3.jar \
    lib/cdk-silent-2.3.jar \
    lib/cdk-smiles-2.3.jar \
    lib/cdk-standard-2.3.jar \
    lib/cdk-valencycheck-2.3.jar \
    lib/commons-net-3.6.jar \
    lib/guava-28.1-jre.jar \
    lib/jdbc-postgresql.jar \
    lib/jnati-core-0.4.jar \
    lib/jnati-deploy-0.4.jar \
    lib/jni-inchi-0.8.jar \
    lib/log4j-1.2.17.jar \
    lib/lucene-core-8.2.0.jar \
    lib/vecmath-1.5.2.jar \
    src/cz/iocb/sachem/load/CompoundLoader.java \
    src/cz/iocb/sachem/load/CompoundUpdater.java \
    src/cz/iocb/sachem/load/ConfigurationProperties.java \
    src/cz/iocb/sachem/load/DrugbankCompoundUpdater.java \
    src/cz/iocb/sachem/load/GeneralCompoundLoader.java \
    src/cz/iocb/sachem/load/PubChemCompoundUpdater.java \
    src/cz/iocb/sachem/lucene/BitSetCollector.java \
    src/cz/iocb/sachem/lucene/FingerprintAnalyzer.java \
    src/cz/iocb/sachem/lucene/FingerprintReader.java \
    src/cz/iocb/sachem/lucene/FingerprintTokenizer.java \
    src/cz/iocb/sachem/lucene/Indexer.java \
    src/cz/iocb/sachem/lucene/Lucene.java \
    src/cz/iocb/sachem/lucene/ScoreHit.java \
    src/cz/iocb/sachem/lucene/Settings.java \
    src/cz/iocb/sachem/lucene/SimilarDocCollector.java \
    src/cz/iocb/sachem/molecule/AromaticityMode.java \
    src/cz/iocb/sachem/molecule/BinaryMoleculeBuilder.java \
    src/cz/iocb/sachem/molecule/BinaryMolecule.java \
    src/cz/iocb/sachem/molecule/BinaryMoleculeSort.java \
    src/cz/iocb/sachem/molecule/MoleculeCreator.java \
    src/cz/iocb/sachem/molecule/Molecule.java \
    src/cz/iocb/sachem/molecule/TautomerMode.java \
    src/cz/iocb/sachem/search/LucyLoader.java \
    src/cz/iocb/sachem/search/QueryFormat.java \
    src/cz/iocb/sachem/search/SimilaritySearch.java \
    src/cz/iocb/sachem/search/SubstructureSearch.java \
    src/cz/iocb/sachem/tautomers/CombinationCountException.java \
    src/cz/iocb/sachem/tautomers/InChIException.java \
    src/cz/iocb/sachem/tautomers/InChI.java \
    src/cz/iocb/sachem/tautomers/InchiTautomerGenerator.java
