create function compound(id in integer) returns varchar language sql as
$$
  select 'http://wifo5-04.informatik.uni-mannheim.de/drugbank/resource/drugs/DB' || right(concat('00000', id), 5);
$$
immutable;


create function compound_inverse(iri in varchar) returns integer language sql as
$$
  select regexp_replace(iri, '^http://wifo5-04.informatik.uni-mannheim.de/drugbank/resource/drugs/DB', '')::integer;
$$
immutable;

--------------------------------------------------------------------------------

create function compound_molfile(id in integer) returns varchar language sql as
$$
  select 'http://wifo5-04.informatik.uni-mannheim.de/drugbank/resource/drugs/DB' || id || '_Molfile';
$$
immutable;


create function compound_molfile_inverse(iri in varchar) returns integer language sql as
$$
  select regexp_replace(iri, '^http://wifo5-04.informatik.uni-mannheim.de/drugbank/resource/drugs/DB([0-9]+)_Molfile$', '\1')::integer;
$$
immutable;
