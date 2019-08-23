create function compound(id in integer) returns varchar language sql as
$$
  select 'http://purl.obolibrary.org/obo/CHEBI_' || id;
$$
immutable parallel safe;


create function compound_inverse(iri in varchar) returns integer language sql as
$$
  select substring(iri, 38)::integer;
$$
immutable parallel safe;

--------------------------------------------------------------------------------

create function compound_molfile(id in integer) returns varchar language sql as
$$
  select 'http://purl.obolibrary.org/obo/CHEBI_' || id || '_Molfile';
$$
immutable parallel safe;


create function compound_molfile_inverse(iri in varchar) returns integer language sql as
$$
  select substring(iri, 38, octet_length(iri) - 45)::integer;
$$
immutable parallel safe;
