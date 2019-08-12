package cz.iocb.sachem.lucene;


enum IndexType
{
    TEXT, POINTS
}


class Settings
{
    static final IndexType indexType = IndexType.TEXT;
    static final boolean useIdTable = true;
    static final boolean useSizeTable = true;
    static final boolean lazyInitialization = true;

    static final String idFieldName = "id";
    static final String subfpFieldName = "subfp";
    static final String simfpFieldName = "simfp";
    static final String simSizeFieldName = "simsz";
}
