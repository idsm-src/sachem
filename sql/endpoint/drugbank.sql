create function compound(id in integer) returns varchar language sql as
$$
  select 'http://wifo5-04.informatik.uni-mannheim.de/drugbank/resource/drugs/DB' || lpad(id::text, 5, '0');
$$
immutable parallel safe;


create function compound_inverse(iri in varchar) returns integer language sql as
$$
  select substring(iri, 70)::integer;
$$
immutable parallel safe;

--------------------------------------------------------------------------------

create function compound_molfile(id in integer) returns varchar language sql as
$$
  select 'http://wifo5-04.informatik.uni-mannheim.de/drugbank/resource/drugs/DB' || lpad(id::text, 5, '0') || '_Molfile';
$$
immutable parallel safe;


create function compound_molfile_inverse(iri in varchar) returns integer language sql as
$$
  select substring(iri, 70, 5)::integer;
$$
immutable parallel safe;
