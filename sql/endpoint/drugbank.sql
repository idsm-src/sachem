create function compound(id in integer) returns varchar language sql as
$$
  select 'http://wifo5-04.informatik.uni-mannheim.de/drugbank/resource/drugs/DB' || right(concat('00000', id), 5);
$$
immutable;
