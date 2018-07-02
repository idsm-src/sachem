create function compound(id in integer) returns varchar language sql as
$$
  select 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID' || id;
$$
immutable;
