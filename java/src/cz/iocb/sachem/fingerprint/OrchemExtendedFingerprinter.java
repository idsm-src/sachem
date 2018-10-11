/*
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work, which includes - but is not limited to -
 * adding the above copyright notice to the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 */
package cz.iocb.sachem.fingerprint;

import java.util.BitSet;
import java.util.Iterator;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import cz.iocb.sachem.fingerprint.bitpos.BitPosApi;
import cz.iocb.sachem.shared.AtomicNumbers;



/**
 * Extends the default {@link OrchemFingerprinter fingerprinter} with further bits of which some are very common. Nice
 * for similarity searching (the more bits the better)
 */
public class OrchemExtendedFingerprinter extends OrchemFingerprinter
{
    public final static int FINGERPRINT_SIZE = BitPosApi.bpExtended.getFingerprintSize();


    @Override
    public int getSize()
    {
        return FINGERPRINT_SIZE;
    }


    /**
     * Returns an extended fingerprint, which is a normal fingerprint with some extra bits at the end.
     *
     * @param molecule to be fingerprinted
     * @return
     */
    @Override
    public BitSet getFingerprint(IAtomContainer molecule) throws CDKException
    {
        return getFingerprint(molecule, 0);
    }


    @Override
    public BitSet getFingerprint(IAtomContainer molecule, int timeout) throws CDKException
    {
        BitSet extendedFingerprint = new BitSet(FINGERPRINT_SIZE);
        Integer bitPos = null;
        String mapKey = null;

        /* Copy the basic fingerprint into the extended fingerprint */
        BitSet basicFingerprint = super.getFingerprint(molecule, timeout);

        for(int i = 0; i < basicFingerprint.size(); i++)
            extendedFingerprint.set(i, basicFingerprint.get(i));

        int atomCount = 0;

        for(IAtom atom : molecule.atoms())
        {
            if(atom.getAtomicNumber() == 0)
                continue;

            /* Fingerprint any element */
            mapKey = atom.getSymbol();
            bitPos = BitPosApi.bpExtended.allElements.get(mapKey);

            if(bitPos != null)
                extendedFingerprint.set(bitPos, true);

            /* Atom count related fingerprinting */
            if(atom.getAtomicNumber() != AtomicNumbers.H)
                atomCount++;

            mapKey = BitPosApi.bpExtended.atomCountMoreThan + atomCount;
            bitPos = BitPosApi.bpExtended.basicStuff.get(mapKey);

            if(bitPos != null)
                extendedFingerprint.set(bitPos, true);
        }

        /* Fingerprint all ring sizes */
        for(IRingSet irs : rslist)
        {
            Iterator<IAtomContainer> ringsetMembers = irs.atomContainers().iterator();

            while(ringsetMembers.hasNext())
            {
                IAtomContainer ring = ringsetMembers.next();
                mapKey = BitPosApi.bp.ringSize + ring.getAtomCount();
                bitPos = BitPosApi.bpExtended.ringSizes.get(mapKey);

                if(bitPos != null)
                    extendedFingerprint.set(bitPos, true);
            }
        }

        /* Fingerprint aromaticity */
        for(IBond bond : molecule.bonds())
        {
            if(bond.getOrder() != null && bond.getFlag(CDKConstants.ISAROMATIC))
            {
                mapKey = BitPosApi.bpExtended.hasAromaticBonds;
                bitPos = BitPosApi.bpExtended.basicStuff.get(mapKey);
                extendedFingerprint.set(bitPos, true);
                break;
            }
        }

        return extendedFingerprint;
    }
}
