CREATE TABLE compounds (
    id                    INT NOT NULL,
    molfile               TEXT NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE orchem_molecules (
    seqid                 INT NOT NULL,
    id                    INT NOT NULL,
    atoms                 BYTEA NOT NULL,
    bonds                 BYTEA NOT NULL,
    PRIMARY KEY (seqid)
);

CREATE TABLE orchem_molecule_counts (
    seqid                 INT NOT NULL,
    molTripleBondCount    SMALLINT NOT NULL,
    molSCount             SMALLINT NOT NULL,
    molOCount             SMALLINT NOT NULL,
    molNCount             SMALLINT NOT NULL,
    molFCount             SMALLINT NOT NULL,
    molClCount            SMALLINT NOT NULL,
    molBrCount            SMALLINT NOT NULL,
    molICount             SMALLINT NOT NULL,
    molCCount             SMALLINT NOT NULL,
    molPCount             SMALLINT NOT NULL,
    PRIMARY KEY (seqid)
);

CREATE TABLE orchem_similarity_fingerprint (
    id                    INT NOT NULL,
    bit_count             SMALLINT NOT NULL,
    fp                    BIGINT[] NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE orchem_substructure_fingerprint (
    seqid                 INT NOT NULL,
    fp                    SMALLINT[] NOT NULL,
    PRIMARY KEY (seqid)
);

CREATE TABLE orchem_substructure_fingerprint_index (
    idx                   SMALLINT NOT NULL,
    bitmap                BYTEA NOT NULL,
    PRIMARY KEY (idx)
);


CREATE INDEX orchem_similarity_fingerprint__bit_count ON orchem_similarity_fingerprint(bit_count);
CREATE INDEX orchem_substructure_fingerprint__fp ON orchem_substructure_fingerprint USING gin(fp);


CREATE FUNCTION "orchem_substructure_search"(varchar, varchar, int, boolean, boolean, boolean, int = 5000) RETURNS SETOF int AS 'libsachem.so','orchem_substructure_search' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "orchem_substructure_gin_search"(varchar, varchar, int, boolean, boolean, boolean, int = 5000) RETURNS SETOF int AS 'libsachem.so','orchem_substructure_gin_search' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "lucy_substructure_search"(varchar, varchar, int, boolean, boolean, boolean, int = 5000) RETURNS SETOF int AS 'libsachem.so','lucy_substructure_search' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "orchem_similarity_search"(varchar, varchar, float4, int) RETURNS TABLE (compound int, score float4) AS 'libsachem.so','orchem_similarity_search' LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION "orchem_substructure_write_indexes"() RETURNS bool AS 'libsachem.so','orchem_substructure_write_indexes' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "lucy_substructure_create_index"() RETURNS bool AS 'libsachem.so','lucy_substructure_create_index' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "lucy_substructure_optimize_index"() RETURNS bool AS 'libsachem.so','lucy_substructure_optimize_index' LANGUAGE C IMMUTABLE STRICT;
