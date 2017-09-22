CREATE TABLE compounds (
    seqid                 INT NOT NULL,
    id                    INT NOT NULL,
    molfile               TEXT NOT NULL,
    PRIMARY KEY (seqid)
);

CREATE TABLE fingerprint_orchem (
    seqid                 INT NOT NULL,
    fp                    BYTEA NOT NULL,
    PRIMARY KEY (seqid)
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


CREATE FUNCTION "orchem_substructure_search"(varchar, varchar, int, boolean, boolean, boolean) RETURNS SETOF int AS 'libpgchem.so','orchem_substructure_search' LANGUAGE C IMMUTABLE STRICT;
