CREATE FUNCTION "orchem_substructure_search"(varchar, varchar, boolean, boolean) RETURNS SETOF int AS 'libpgchem.so','orchem_substructure_search' LANGUAGE C IMMUTABLE STRICT;
