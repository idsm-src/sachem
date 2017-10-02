CREATE TABLE compounds (
    id                    INT NOT NULL,
    molfile               TEXT NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE fingerprint_orchem (
    id                    INT NOT NULL,
    bit_count             SMALLINT NOT NULL,
    fp                    BIGINT[] NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE molecule_counts (
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

CREATE TABLE molecules (
    seqid                 INT NOT NULL,
    id                    INT NOT NULL,
    atoms                 BYTEA NOT NULL,
    bonds                 BYTEA NOT NULL,
    PRIMARY KEY (seqid)
);

CREATE TABLE fingerprint_orchem_index (
    idx                   SMALLINT NOT NULL,
    bitmap                BYTEA NOT NULL,
    PRIMARY KEY (idx)
);


CREATE INDEX fingerprint_orchem__bit_count ON fingerprint_orchem(bit_count);


CREATE FUNCTION "orchem_substructure_search"(varchar, varchar, int, boolean, boolean, boolean) RETURNS SETOF int AS 'libpgchem.so','orchem_substructure_search' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "orchem_similarity_search"(varchar, varchar, float4, int) RETURNS TABLE (compound int, score float4) AS 'libpgchem.so','orchem_similarity_search' LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION "orchem_substructure_write_indexes"() RETURNS bool AS 'libpgchem.so','orchem_substructure_write_indexes' LANGUAGE C IMMUTABLE STRICT;
