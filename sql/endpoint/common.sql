create function queryFormat_inverse(iri in varchar) returns integer language plpgsql as
$$
  declare value record;
  begin
    select into value 
      case iri
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/SMILES' then 0
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/MolFile' then 1
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/RGroup' then 2
      end::integer AS retval;
    return value.retval;
  end;
$$
immutable;


create function graphMode_inverse(iri in varchar) returns integer language plpgsql as
$$
  declare value record;
  begin
    select into value 
      case iri
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/substructureSearch' then 0
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/exactSearch' then 1
      end::integer AS retval;
    return value.retval;
  end;
$$
immutable;


create function chargeMode_inverse(iri in varchar) returns integer language plpgsql as
$$
  declare value record;
  begin
    select into value 
      case iri
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/ignoreCharges' then 0
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/defaultChargeAsZero' then 1
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/defaultChargeAsAny' then 2
      end::integer AS retval;
    return value.retval;
  end;
$$
immutable;


create function isotopeMode_inverse(iri in varchar) returns integer language plpgsql as
$$
  declare value record;
  begin
    select into value 
      case iri
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/ignoreIsotopes' then 0
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/defaultIsotopeAsStandard' then 1
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/defaultIsotopeAsAny' then 2
      end::integer AS retval;
    return value.retval;
  end;
$$
immutable;


create function stereoMode_inverse(iri in varchar) returns integer language plpgsql as
$$
  declare value record;
  begin
    select into value 
      case iri
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/ignoreStrereo' then 0
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/strictStereo' then 1
      end::integer AS retval;
    return value.retval;
  end;
$$
immutable;


create function tautomerMode_inverse(iri in varchar) returns integer language plpgsql as
$$
  declare value record;
  begin
    select into value 
      case iri
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/ignoreTautomers' then 0
        when 'http://bioinfo.uochb.cas.cz/sparql-endpoint/sachem/inchiTautomers' then 1
      end::integer AS retval;
    return value.retval;
  end;
$$
immutable;
