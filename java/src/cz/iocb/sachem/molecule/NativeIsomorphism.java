package cz.iocb.sachem.molecule;

import java.nio.ByteBuffer;



public class NativeIsomorphism
{
    @SuppressWarnings("serial")
    public static class QueryCancelException extends RuntimeException
    {
    }

    @SuppressWarnings("serial")
    public static class IterationLimitExceededException extends Exception
    {
    }


    private final ByteBuffer implementation;


    public NativeIsomorphism(byte[] query, boolean[] restH, SearchMode graphMode, ChargeMode chargeMode,
            IsotopeMode isotopeMode, StereoMode stereoMode)
    {
        implementation = create(query, restH, graphMode.ordinal(), chargeMode.ordinal(), isotopeMode.ordinal(),
                stereoMode.ordinal());
    }


    public float match(byte[] query, int limit) throws IterationLimitExceededException
    {
        return match(implementation, query, limit);
    }


    private static native float match(ByteBuffer implementation, byte[] query, int limit)
            throws IterationLimitExceededException;


    private static native ByteBuffer create(byte[] query, boolean[] restH, int graphMode, int chargeMode,
            int isotopeMode, int stereoMode);
}
