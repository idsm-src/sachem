-- complain if script is sourced in psql, rather than via CREATE EXTENSION
\echo Use "CREATE EXTENSION sachem" to load this file. \quit


CREATE TYPE search_mode AS ENUM ('SUBSTRUCTURE', 'EXACT');
CREATE TYPE charge_mode AS ENUM ('IGNORE', 'DEFAULT_AS_UNCHARGED', 'DEFAULT_AS_ANY');
CREATE TYPE isotope_mode AS ENUM ('IGNORE', 'DEFAULT_AS_STANDARD', 'DEFAULT_AS_ANY');
CREATE TYPE stereo_mode AS ENUM ('IGNORE', 'STRICT');
CREATE TYPE aromaticity_mode AS ENUM ('PRESERVE', 'DETECT', 'AUTO');
CREATE TYPE tautomer_mode AS ENUM ('IGNORE', 'INCHI');
CREATE TYPE query_format AS ENUM ('UNSPECIFIED', 'SMILES', 'MOLFILE', 'RGROUP');


CREATE TABLE configuration (
    id              SERIAL NOT NULL,
    index_name      VARCHAR NOT NULL UNIQUE CHECK (index_name ~ '^[a-zA-Z0-9_]+$'),
    schema_name     VARCHAR NOT NULL CHECK (char_length(schema_name) > 0),
    table_name      VARCHAR NOT NULL CHECK (char_length(table_name) > 0),
    id_column       VARCHAR NOT NULL CHECK (char_length(id_column) > 0),
    molfile_column  VARCHAR NOT NULL CHECK (char_length(molfile_column) > 0),
    threads         INT NOT NULL CHECK (char_length(molfile_column) > 0),
    segments        INT NOT NULL CHECK (char_length(molfile_column) > 0),
    buffered_docs   INT NOT NULL CHECK (char_length(molfile_column) >= 0),
    buffer_size     FLOAT4 NOT NULL CHECK (char_length(molfile_column) >= 0),
    version         INT NOT NULL,
    PRIMARY KEY (id)
);


CREATE TABLE compound_audit (
    index                 INT NOT NULL REFERENCES configuration(id),
    id                    INT NOT NULL,
    stored                BOOLEAN NOT NULL,
    PRIMARY KEY (index, id)
);


CREATE TABLE compound_errors (
    id                    SERIAL NOT NULL,
    timestamp             TIMESTAMPTZ NOT NULL DEFAULT now(),
    index                 INT NOT NULL REFERENCES configuration(id),
    compound              INT NOT NULL,
    message               TEXT NOT NULL,
    PRIMARY KEY (id)
);


CREATE TABLE compound_sources (
    id                    SERIAL NOT NULL,
    index                 INT NOT NULL REFERENCES configuration(id),
    name                  TEXT NOT NULL,
    size                  BIGINT NOT NULL,
    timestamp             TIMESTAMPTZ,
    PRIMARY KEY (id)
);


CREATE TABLE compound_stats (
    index                 INT NOT NULL REFERENCES configuration(id),
    version               TEXT,
    checkdate             TIMESTAMPTZ NOT NULL,
    PRIMARY KEY (index)
);


GRANT USAGE ON SEQUENCE compound_sources_id_seq TO PUBLIC;

GRANT SELECT ON TABLE configuration TO PUBLIC;
GRANT SELECT ON TABLE compound_sources TO PUBLIC;
GRANT INSERT ON TABLE compound_sources TO PUBLIC;
GRANT UPDATE ON TABLE compound_sources TO PUBLIC;
GRANT DELETE ON TABLE compound_sources TO PUBLIC;
GRANT SELECT ON TABLE compound_stats TO PUBLIC;
GRANT INSERT ON TABLE compound_stats TO PUBLIC;
GRANT UPDATE ON TABLE compound_stats TO PUBLIC;
GRANT SELECT ON TABLE compound_errors TO PUBLIC;


CREATE FUNCTION "index_size"(varchar) RETURNS int AS 'MODULE_PATHNAME' LANGUAGE C IMMUTABLE;
CREATE FUNCTION "substructure_search"(varchar, varchar, search_mode = 'SUBSTRUCTURE', charge_mode = 'DEFAULT_AS_ANY', isotope_mode = 'IGNORE', stereo_mode = 'IGNORE', aromaticity_mode = 'AUTO', tautomer_mode = 'IGNORE', query_format = 'UNSPECIFIED', int = -1, boolean = false, int = 0) RETURNS TABLE (compound int, score float4) AS 'MODULE_PATHNAME' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "similarity_search"(varchar, varchar, float4 = 0.85, int = 1, aromaticity_mode = 'AUTO', tautomer_mode = 'IGNORE', query_format = 'UNSPECIFIED', int = -1, boolean = false) RETURNS TABLE (compound int, score float4) AS 'MODULE_PATHNAME' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "sync_data"(varchar, boolean = false, boolean = true) RETURNS void AS 'MODULE_PATHNAME' LANGUAGE C IMMUTABLE STRICT SECURITY DEFINER;
CREATE FUNCTION "cleanup"(varchar) RETURNS void AS 'MODULE_PATHNAME' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION "segments"(varchar) RETURNS TABLE (name varchar, molecules int, deletes int, size bigint) AS 'MODULE_PATHNAME' LANGUAGE C IMMUTABLE STRICT;


CREATE FUNCTION "add_index"(index_name varchar, schema_name varchar, table_name varchar, id_column varchar = 'id', molfile_column varchar = 'molfile', threads int = 4, segments int = 4, buffered_docs int = 1000, buffer_size float8 = 64) RETURNS void AS $$
DECLARE
    idx int;
BEGIN
	INSERT INTO sachem.configuration (index_name, schema_name, table_name, id_column, molfile_column, threads, segments, buffered_docs, buffer_size, version) VALUES (index_name, schema_name, table_name, id_column, molfile_column, threads, segments, buffered_docs, buffer_size, 0);

	SELECT id INTO idx FROM sachem.configuration AS tbl WHERE tbl.index_name = "add_index".index_name;
	
	schema_name := quote_ident(schema_name);
	table_name := quote_ident(table_name);
    id_column := quote_ident(id_column);
    molfile_column := quote_ident(molfile_column);	
    
	EXECUTE 'CREATE FUNCTION sachem."' || index_name || '_compound_audit"() RETURNS TRIGGER AS
	$body$
	BEGIN
	  IF (TG_OP = ''INSERT'') THEN
	    INSERT INTO sachem.compound_audit (index, id, stored) VALUES (' || idx || ', NEW.' || id_column || ', true)
	        ON CONFLICT (index, id) DO UPDATE SET index=EXCLUDED.index, id=EXCLUDED.id, stored=true;
	    RETURN NEW;
	  ELSIF (TG_OP = ''UPDATE'') THEN
	    IF (OLD.' || molfile_column || ' != NEW.' || molfile_column || ') THEN
	        INSERT INTO sachem.compound_audit (index, id, stored) VALUES (' || idx || ', NEW.' || id_column || ', true)
	            ON CONFLICT (index, id) DO UPDATE SET index=EXCLUDED.index, id=EXCLUDED.id, stored=true;
	    END IF;
	    RETURN NEW;
	  ELSIF (TG_OP = ''DELETE'') THEN
	    INSERT INTO sachem.compound_audit (index, id, stored) VALUES (' || idx || ', OLD.' || id_column || ', false)
	        ON CONFLICT (index, id) DO UPDATE SET index=EXCLUDED.index, id=EXCLUDED.id, stored=false;
	    RETURN OLD;
	  ELSIF (TG_OP = ''TRUNCATE'') THEN
	    INSERT INTO sachem.compound_audit SELECT ' || idx || ', ' || id_column || ', false FROM ' || schema_name || '.' || table_name || '
	        ON CONFLICT (index, id) DO UPDATE SET index=EXCLUDED.index, id=EXCLUDED.id, stored=false;
	    RETURN NULL;
	  END IF;
	END;
	$body$ LANGUAGE plpgsql SECURITY DEFINER';
	
    EXECUTE 'CREATE TRIGGER "' || index_name || '_sachem_compound_audit" AFTER INSERT OR UPDATE OR DELETE ON ' ||
        schema_name || '.' || table_name || ' FOR EACH ROW EXECUTE PROCEDURE sachem."' || index_name || '_compound_audit"()';
	
	EXECUTE 'CREATE TRIGGER "' || index_name || '_sachem_truncate_compound_audit" BEFORE TRUNCATE ON ' ||
	    schema_name || '.' || table_name ||  ' FOR EACH STATEMENT EXECUTE PROCEDURE sachem."' || index_name || '_compound_audit"()';
    
	EXECUTE 'INSERT INTO sachem.compound_audit (index, id, stored) SELECT ' || idx || ', ' || id_column || ', true FROM ' || schema_name || '.' || table_name;
END
$$ LANGUAGE PLPGSQL STRICT SECURITY DEFINER;


CREATE FUNCTION "remove_index"(index_name varchar) RETURNS void AS $$
DECLARE
    idx int;
    schema_name varchar;
    table_name varchar;
BEGIN
	SELECT id INTO idx FROM sachem.configuration AS tbl WHERE tbl.index_name = "remove_index".index_name;
    SELECT quote_ident(tbl.schema_name) INTO schema_name FROM sachem.configuration AS tbl WHERE tbl.index_name = "remove_index".index_name;
    SELECT quote_ident(tbl.table_name) INTO table_name FROM sachem.configuration AS tbl WHERE tbl.index_name = "remove_index".index_name;
	
	DELETE FROM sachem.compound_audit WHERE index = idx;
	DELETE FROM sachem.compound_errors WHERE index = idx;
	DELETE FROM sachem.compound_sources WHERE index = idx;
	DELETE FROM sachem.compound_stats WHERE index = idx;
	DELETE FROM sachem.configuration AS tbl WHERE tbl.id = idx;
	
	EXECUTE 'DROP TRIGGER "' || index_name || '_sachem_compound_audit" ON ' || schema_name || '.' || table_name;
    EXECUTE 'DROP TRIGGER "' || index_name || '_sachem_truncate_compound_audit" ON ' || schema_name || '.' || table_name;
    EXECUTE 'DROP FUNCTION sachem."' || index_name || '_compound_audit"()';
    
    PERFORM sachem.cleanup(index_name);
END
$$ LANGUAGE PLPGSQL STRICT SECURITY DEFINER;
