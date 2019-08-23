create function compound(id in integer) returns varchar language sql as
$$
  select 'http://rdf.ebi.ac.uk/resource/chembl/molecule/CHEMBL' || id;
$$
immutable parallel safe;


create function compound_inverse(iri in varchar) returns integer language sql as
$$
  select substring(iri, 53)::integer;
$$
immutable parallel safe;

--------------------------------------------------------------------------------

create function compound_molfile(id in integer) returns varchar language sql as
$$
  select 'http://rdf.ebi.ac.uk/resource/chembl/molecule/CHEMBL' || id || '_Molfile';
$$
immutable parallel safe;


create function compound_molfile_inverse(iri in varchar) returns integer language sql as
$$
  select substring(iri, 53, octet_length(iri) - 60)::integer;
$$
immutable parallel safe;
