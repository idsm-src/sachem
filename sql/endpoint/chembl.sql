create function compound(id in integer) returns varchar language sql as
$$
  select 'http://rdf.ebi.ac.uk/resource/chembl/molecule/CHEMBL' || id;
$$
immutable;


create function compound_inverse(iri in varchar) returns integer language sql as
$$
  select regexp_replace(iri, '^http://rdf.ebi.ac.uk/resource/chembl/molecule/CHEMBL', '')::integer;
$$
immutable;

--------------------------------------------------------------------------------

create function compound_molfile(id in integer) returns varchar language sql as
$$
  select 'http://rdf.ebi.ac.uk/resource/chembl/molecule/CHEMBL' || id || '_Molfile';
$$
immutable;


create function compound_molfile_inverse(iri in varchar) returns integer language sql as
$$
  select regexp_replace(iri, '^http://rdf.ebi.ac.uk/resource/chembl/molecule/CHEMBL([0-9]+)_Molfile$', '\1')::integer;
$$
immutable;
