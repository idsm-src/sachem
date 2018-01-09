create function chembl(id in integer) returns varchar language sql as
$$
  select 'http://rdf.ebi.ac.uk/resource/chembl/molecule/CHEMBL' || id;
$$
immutable;
