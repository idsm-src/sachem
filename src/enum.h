#ifndef ENUM_H_
#define ENUM_H_

#include <postgres.h>
#include <catalog/pg_enum.h>
#include <utils/syscache.h>
#include <access/htup_details.h>
#include "java.h"


typedef struct
{
    Oid oid;
    jobject object;
}
EnumValue;


static inline jobject ConvertEnumValue(EnumValue *table, Oid oid)
{
    for(int i = 0; table[i].oid != InvalidOid; i++)
        if(table[i].oid == oid)
            return table[i].object;

    return NULL;
}


static inline Oid LookupExplicitEnumType(const Oid spaceid, const char *name)
{
    #if PG_VERSION_NUM >= 120000
    Oid typoid = GetSysCacheOid2(TYPENAMENSP, Anum_pg_type_oid, PointerGetDatum(name), ObjectIdGetDatum(spaceid));
    #else
    Oid typoid = GetSysCacheOid2(TYPENAMENSP, PointerGetDatum(name), ObjectIdGetDatum(spaceid));
    #endif

    if(!OidIsValid(typoid))
        elog(ERROR, "cannot find enum type '%s'", name);

    return typoid;
}


static inline Oid LookupExplicitEnumValue(const Oid typoid, const char *name)
{
    HeapTuple tup = SearchSysCache2(ENUMTYPOIDNAME, ObjectIdGetDatum(typoid), CStringGetDatum(name));

    if(!HeapTupleIsValid(tup))
        elog(ERROR, "cannot find enum value '%s'", name);

    #if PG_VERSION_NUM >= 120000
    Oid valueid = ((Form_pg_enum) GETSTRUCT(tup))->oid;
    #else
    Oid valueid = HeapTupleGetOid(tup);
    #endif

    ReleaseSysCache(tup);

    return valueid;
}


static inline jobject LookupJavaEnumValue(jclass enumClass, const char *name, const char *sig)
{
    jfieldID fieldId = (*env)->GetStaticFieldID(env, enumClass, name, sig);
    java_check_exception(__func__);

    jobject fieldValue = (jobject) (*env)->NewGlobalRef(env, (*env)->GetStaticObjectField(env, enumClass, fieldId));
    java_check_exception(__func__);

    return fieldValue;
}

#endif /* ENUM_H_ */
