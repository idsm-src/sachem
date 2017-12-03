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
    id                    INT NOT NULL,
    counts                SMALLINT[] NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE orchem_fingerprint (
    id                    INT NOT NULL,
    bit_count             SMALLINT NOT NULL,
    fp                    BIGINT[] NOT NULL,
    PRIMARY KEY (id)
);


CREATE INDEX orchem_fingerprint__bit_count ON orchem_fingerprint(bit_count);


CREATE FUNCTION "orchem_substructure_search"(varchar, varchar, int, boolean, boolean, boolean, int = 5000) RETURNS SETOF int AS 'libsachem.so','orchem_substructure_search' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "lucy_substructure_search"(varchar, varchar, int, boolean, boolean, boolean, int = 5000) RETURNS SETOF int AS 'libsachem.so','lucy_substructure_search' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "orchem_similarity_search"(varchar, varchar, float4, int) RETURNS TABLE (compound int, score float4) AS 'libsachem.so','orchem_similarity_search' LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION "orchem_substructure_write_indexes"() RETURNS bool AS 'libsachem.so','orchem_substructure_write_indexes' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "orchem_load_data"() RETURNS int AS 'libsachem.so','orchem_load_data' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "lucy_substructure_create_index"() RETURNS bool AS 'libsachem.so','lucy_substructure_create_index' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "lucy_substructure_optimize_index"() RETURNS bool AS 'libsachem.so','lucy_substructure_optimize_index' LANGUAGE C IMMUTABLE STRICT;
