create function compound(id in integer) returns varchar language sql as
$$
  select 'http://purl.obolibrary.org/obo/CHEBI_' || id;
$$
immutable;


create function compound_inverse(iri in varchar) returns integer language sql as
$$
  select regexp_replace(iri, '^http://purl.obolibrary.org/obo/CHEBI_', '')::integer;
$$
immutable;

--------------------------------------------------------------------------------

create function compound_molfile(id in integer) returns varchar language sql as
$$
  select 'http://purl.obolibrary.org/obo/CHEBI_' || id || '_Molfile';
$$
immutable;


create function compound_molfile_inverse(iri in varchar) returns integer language sql as
$$
  select regexp_replace(iri, '^http://purl.obolibrary.org/obo/CHEBI_([0-9]+)_Molfile$', '\1')::integer;
$$
immutable;
