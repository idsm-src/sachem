-- complain if script is sourced in psql, rather than via CREATE EXTENSION
\echo Use "CREATE EXTENSION sachem_orchem" to load this file. \quit


CREATE TABLE compounds (
    id                    INT NOT NULL,
    molfile               TEXT NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE compound_sources (
    id                    SERIAL NOT NULL,
    name                  TEXT NOT NULL,
    size                  BIGINT NOT NULL,
    timestamp             TIMESTAMPTZ,
    PRIMARY KEY (id)
);

CREATE TABLE compound_stats (
    id                    SMALLINT,
    version               TEXT,
    size                  BIGINT NOT NULL,
    checkdate             TIMESTAMPTZ NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE orchem_compound_audit (
    id                    INT NOT NULL,
    stored                BOOLEAN NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE orchem_index (
    id                    INT NOT NULL,
    path                  TEXT NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE orchem_molecules (
    id                    INT NOT NULL,
    seqid                 INT NOT NULL,
    atoms                 BYTEA NOT NULL,
    bonds                 BYTEA NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE orchem_molecule_errors (
    id                    SERIAL NOT NULL,
    timestamp             TIMESTAMPTZ NOT NULL DEFAULT now(),
    compound              INT NOT NULL,
    message               TEXT NOT NULL,
    PRIMARY KEY (id)
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


CREATE INDEX orchem_molecules__seqid ON orchem_molecules(seqid);
CREATE INDEX orchem_fingerprint__bit_count ON orchem_fingerprint(bit_count);


GRANT SELECT ON TABLE compounds TO PUBLIC;
GRANT INSERT ON TABLE compounds TO PUBLIC;
GRANT UPDATE ON TABLE compounds TO PUBLIC;
GRANT DELETE ON TABLE compounds TO PUBLIC;
GRANT TRUNCATE ON TABLE compounds TO PUBLIC;
GRANT SELECT ON TABLE compound_sources TO PUBLIC;
GRANT INSERT ON TABLE compound_sources TO PUBLIC;
GRANT UPDATE ON TABLE compound_sources TO PUBLIC;
GRANT DELETE ON TABLE compound_sources TO PUBLIC;
GRANT TRUNCATE ON TABLE compound_sources TO PUBLIC;
GRANT SELECT ON SEQUENCE compound_sources_id_seq TO PUBLIC;
GRANT USAGE ON SEQUENCE compound_sources_id_seq TO PUBLIC;
GRANT SELECT ON TABLE compound_stats TO PUBLIC;
GRANT INSERT ON TABLE compound_stats TO PUBLIC;
GRANT UPDATE ON TABLE compound_stats TO PUBLIC;
GRANT DELETE ON TABLE compound_stats TO PUBLIC;
GRANT TRUNCATE ON TABLE compound_stats TO PUBLIC;
GRANT SELECT ON TABLE orchem_molecule_errors TO PUBLIC;


CREATE FUNCTION "orchem_substructure_search"(varchar, varchar, int, boolean, boolean, boolean, int = 5000) RETURNS SETOF int AS 'MODULE_PATHNAME','orchem_substructure_search' LANGUAGE C IMMUTABLE STRICT SECURITY DEFINER;
CREATE FUNCTION "orchem_similarity_search"(varchar, varchar, float4, int) RETURNS TABLE (compound int, score float4) AS 'MODULE_PATHNAME','orchem_similarity_search' LANGUAGE C IMMUTABLE STRICT SECURITY DEFINER;
CREATE FUNCTION "orchem_sync_data"(boolean = false) RETURNS void AS 'MODULE_PATHNAME','orchem_sync_data' LANGUAGE C IMMUTABLE STRICT SECURITY DEFINER;


CREATE FUNCTION orchem_compound_audit() RETURNS TRIGGER AS
$body$
BEGIN
  IF (TG_OP = 'INSERT') THEN
    INSERT INTO orchem_compound_audit (id, stored) VALUES (NEW.id, true)
        ON CONFLICT (id) DO UPDATE SET id=EXCLUDED.id, stored=true;
    RETURN NEW;
  ELSIF (TG_OP = 'UPDATE') THEN
    IF (OLD.molfile != NEW.molfile) THEN
        INSERT INTO orchem_compound_audit (id, stored) VALUES (NEW.id, true)
            ON CONFLICT (id) DO UPDATE SET id=EXCLUDED.id, stored=true;
    END IF;
    RETURN NEW;
  ELSIF (TG_OP = 'DELETE') THEN
    INSERT INTO orchem_compound_audit (id, stored) VALUES (OLD.id, false)
        ON CONFLICT (id) DO UPDATE SET id=EXCLUDED.id, stored=false;
    RETURN OLD;
  ELSIF (TG_OP = 'TRUNCATE') THEN
    INSERT INTO orchem_compound_audit SELECT id, false FROM compounds
        ON CONFLICT (id) DO UPDATE SET id=EXCLUDED.id, stored=false;
    RETURN NULL;
  END IF;
END;
$body$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE TRIGGER orchem_compound_audit AFTER INSERT OR UPDATE OR DELETE ON compounds
    FOR EACH ROW EXECUTE PROCEDURE orchem_compound_audit();

CREATE TRIGGER orchem_truncate_compounds_audit BEFORE TRUNCATE ON compounds
    FOR EACH STATEMENT EXECUTE PROCEDURE orchem_compound_audit();
