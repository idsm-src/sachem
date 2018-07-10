create table class_bases
(
    id       integer,
    iri      varchar not null,
    primary key(id),
    unique(iri)
);


create table property_bases
(
    id       integer,
    iri      varchar not null,
    primary key(id),
    unique(iri)
);

--------------------------------------------------------------------------------

insert into class_bases(id, iri) values (0, 'http://semanticscience.org/resource/SIO_011120');
insert into property_bases(id, iri) values (0, 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type');
insert into property_bases(id, iri) values (1, 'http://semanticscience.org/resource/is-attribute-of');
insert into property_bases(id, iri) values (2, 'http://semanticscience.org/resource/has-value');

--------------------------------------------------------------------------------

create function class(id in integer) returns varchar language sql as
$$
  select iri from class_bases where class_bases.id = class.id;
$$
immutable;


create function class_inverse(iri in varchar) returns integer language sql as
$$
  select id from class_bases where class_bases.iri = class_inverse.iri;
$$
immutable;


create function property(id in integer) returns varchar language sql as
$$
  select iri from property_bases where property_bases.id = property.id;
$$
immutable;


create function property_inverse(iri in varchar) returns integer language sql as
$$
  select id from property_bases where property_bases.iri = property_inverse.iri;
$$
immutable;
