create function chebi(id in integer) returns varchar language sql as
$$
  select 'http://purl.obolibrary.org/obo/CHEBI_' || id;
$$
immutable;
