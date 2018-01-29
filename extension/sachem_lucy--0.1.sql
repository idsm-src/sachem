-- complain if script is sourced in psql, rather than via CREATE EXTENSION
\echo Use "CREATE EXTENSION sachem_lucy" to load this file. \quit


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

CREATE TABLE sachem_compound_audit (
    id                    INT NOT NULL,
    stored                BOOLEAN NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE sachem_molecules (
    id                    INT NOT NULL,
    molecule              BYTEA NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE sachem_molecule_errors (
    id                    SERIAL NOT NULL,
    timestamp             TIMESTAMPTZ NOT NULL DEFAULT now(),
    compound              INT NOT NULL,
    message               TEXT NOT NULL,
    PRIMARY KEY (id)
);


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
GRANT SELECT ON TABLE sachem_molecule_errors TO PUBLIC;


CREATE FUNCTION "sachem_substructure_search"(varchar, int, int = 0, int = 0, int = 2, int = 0, int = 0, int = 0, int = 5000) RETURNS SETOF int AS 'MODULE_PATHNAME','lucy_substructure_search' LANGUAGE C IMMUTABLE STRICT SECURITY DEFINER;
CREATE FUNCTION "sachem_sync_data"(boolean = false) RETURNS void AS 'MODULE_PATHNAME','lucy_sync_data' LANGUAGE C IMMUTABLE STRICT SECURITY DEFINER;


CREATE FUNCTION sachem_compound_audit() RETURNS TRIGGER AS
$body$
BEGIN
  IF (TG_OP = 'INSERT') THEN
    INSERT INTO sachem_compound_audit (id, stored) VALUES (NEW.id, true)
        ON CONFLICT (id) DO UPDATE SET id=EXCLUDED.id, stored=true;
    RETURN NEW;
  ELSIF (TG_OP = 'UPDATE') THEN
    IF (OLD.molfile != NEW.molfile) THEN
        INSERT INTO sachem_compound_audit (id, stored) VALUES (NEW.id, true)
            ON CONFLICT (id) DO UPDATE SET id=EXCLUDED.id, stored=true;
    END IF;
    RETURN NEW;
  ELSIF (TG_OP = 'DELETE') THEN
    INSERT INTO sachem_compound_audit (id, stored) VALUES (OLD.id, false)
        ON CONFLICT (id) DO UPDATE SET id=EXCLUDED.id, stored=false;
    RETURN OLD;
  ELSIF (TG_OP = 'TRUNCATE') THEN
    INSERT INTO sachem_compound_audit SELECT id, false FROM compounds
        ON CONFLICT (id) DO UPDATE SET id=EXCLUDED.id, stored=false;
    RETURN NULL;
  END IF;
END;
$body$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE TRIGGER sachem_compound_audit AFTER INSERT OR UPDATE OR DELETE ON compounds
    FOR EACH ROW EXECUTE PROCEDURE sachem_compound_audit();

CREATE TRIGGER sachem_truncate_compounds_audit BEFORE TRUNCATE ON compounds
    FOR EACH STATEMENT EXECUTE PROCEDURE sachem_compound_audit();
