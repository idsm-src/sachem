/*
 * Copyright (C) 2017-2017 Jakub Galgonek   galgonek@uochb.cas.cz
 * Copyright (C) 2008-2009 Mark Rijnbeek    markr@ebi.ac.uk
 *
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
package cz.iocb.sachem.fingerprint.bitpos;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IBond;



/**
 * Holds maps with fingerprint information used to determine which bit number to set for which chemical aspect. Access
 * to this class should be done through singleton {@link BitPosApi}
 *
 * Several maps are set up, each containing a pair (information,bit position). For example, the map elementCntBits holds
 * data related to element counts. It may have an entry (C32,n) with n being some number - this data will be used by
 * {@link cz.iocb.sachem.fingerprint.OrchemFingerprinter} to determine which bit position (n) to set when a compound has
 * 32 or more Carbon atoms. And so on.. some maps are more elaborate than others but should have reasonably
 * self-explaining content.
 */
public class BitPositions
{
    /* Map for fingerprinting element counts */
    public final Map<String, Integer> elemCntBits = new HashMap<String, Integer>();
    public final Map<String, String> elemFingerprinted = new HashMap<String, String>();

    /* Map for atom pair fingerprinting */
    public final Map<String, Integer> atomPairBits = new HashMap<String, Integer>();

    /* Map for element neighbour fingerprinting */
    public final Map<Integer, List<Neighbour>> neighbourBits = new HashMap<Integer, List<Neighbour>>();

    /* Map for smarts pattern fingerprinting */
    public final Map<String, Integer> smilesPatternBits = new HashMap<String, Integer>();

    /* Map for fingerprinting ringsets with multiple (!) rings */
    public final Map<String, Integer> ringSetBits = new HashMap<String, Integer>();

    /* Map for fingerprinting ring aspects */
    public final Map<String, Integer> ringBits = new HashMap<String, Integer>();

    /* Map for fingerprinting ring sizes */
    public final Map<String, Integer> ringSizes = new HashMap<String, Integer>();

    /* Map for fingerprinting carbon trails */
    public final Map<String, Integer> carbonTrails = new HashMap<String, Integer>();
    public final int MAX_CC_TRAIL_DEPTH = 20;

    /* Map for fingerprinting ring layout */
    public final Map<String, Integer> ringLayout = new HashMap<String, Integer>();


    /* Reserved prefix strings for readable labeling of bit position dump */
    public final String ringsetCountTotalPrefix = new String("RsTot");
    public final String ringsetCountPerSet = new String("RsSet");
    public final String ringsetPairPrefix = new String("RsPair");
    public final String ringsetCountConnectedPrefix = new String("RsNeigh");
    public final String ringsetCountHexPrefix = new String("RsHex");
    public final String ringsetCountPentPrefix = new String("RsPent");
    public final String ringPrefixRing = new String("ring_");
    public final String ringPrefixAny = new String("any_");
    public final String ringPrefixArom = new String("arom_");
    public final String ringPrefixNonArom = new String("nonArom_");
    public final String ringPrefixNitro = new String("N_");
    public final String ringPrefixHetero = new String("Het_");
    public final String ringPrefixCarbonOnly = new String("CarbOnly_");
    public final String cTrailPrefix = new String("c-c_Trail");
    public final String hexRingAromNeighPrefix = new String("hexRingAromNeigh_");
    public final String hexRingNonAromNeighPrefix = new String("hexRingNonAromNeigh_");
    public final String hexRingMixNeighPrefix = new String("hexRingMixAromNeigh_");
    public final String ringNonaromDoubleBond = new String("ringDoubleBondNonArom_");
    public final String otherElem = "OtherElement";
    public final String ringSize = new String("ringSize_");

    public final int HASH_OFFSET = 95;

    int bitpos = 0;


    /**
     * Returns the size of the fingerprint = the value of bitPos as incremented in the constructor
     */
    public int getFingerprintSize()
    {
        return bitpos + 1;
    }


    /**
     * Constructor sets up the maps with fingerprinting data.
     */
    public BitPositions()
    {
        /*
         * Position 0 is reserved !
         * bit0 is ALWAYS set and thus prevents division by zero in the similarity search (weakness of algorithm)
         *
         */

        bitpos = HASH_OFFSET;


        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("C", IBond.Order.SINGLE, false) }));
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + "2", ++bitpos);
        elemCntBits.put("O3", ++bitpos);
        elemCntBits.put("C20", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixAny + "1", ++bitpos);
        ringSetBits.put(ringsetCountTotalPrefix + "1", ++bitpos);
        elemCntBits.put("N3", ++bitpos);
        smilesPatternBits.put("C-C-N-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("N", null, false), new Neighbour("O", null, false) }));
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixCarbonOnly + "2", ++bitpos);
        smilesPatternBits.put("C:C-C-C", ++bitpos);
        smilesPatternBits.put("C-C-C-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("N", IBond.Order.SINGLE, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("O", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, false) }));
        ringBits.put("6" + ringPrefixRing + ringPrefixAny + "3", ++bitpos);
        smilesPatternBits.put("N-C-C-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, false), new Neighbour("C", null, false) }));
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                        new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("N", IBond.Order.SINGLE, false),
                        new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("O=C-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, false), new Neighbour("C", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, false), new Neighbour("N", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("N", null, true) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, false),
                new Neighbour("C", null, true), new Neighbour("C", null, true) }));
        smilesPatternBits.put("O=C-N-C-C", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + "1", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("N", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("O", null, false), new Neighbour("C", null, true) }));
        elemCntBits.put("C24", ++bitpos);
        elemCntBits.put("N4", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, false), new Neighbour("O", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("O", null, false), new Neighbour("O", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("O", null, false), new Neighbour("C", null, true), new Neighbour("C", null, true) }));
        smilesPatternBits.put("C-C-C-C:C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + "1", ++bitpos);
        smilesPatternBits.put("O=C-C-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("O", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("O-C-C-C-C", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "1", ++bitpos);
        smilesPatternBits.put("N-C-C:C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("C", null, true), new Neighbour("N", null, false) }));
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + "3", ++bitpos);
        ringSetBits.put(ringsetPairPrefix + "5_6", ++bitpos);
        smilesPatternBits.put("C-C-C-C-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("O", null, false), new Neighbour("O", null, false) }));
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "1", ++bitpos);
        elemCntBits.put("O5", ++bitpos);
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                        new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("O", IBond.Order.SINGLE, false),
                        new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("C-N-C:C", ++bitpos);
        smilesPatternBits.put("C:C-O-C", ++bitpos);
        smilesPatternBits.put("C-N-C-C:C", ++bitpos);
        smilesPatternBits.put("N-C-C-N", ++bitpos);
        elemCntBits.put("S1", ++bitpos);
        smilesPatternBits.put("N-C-C-N-C", ++bitpos);
        smilesPatternBits.put("O=C-C:C", ++bitpos);
        smilesPatternBits.put("O-C-C-N", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("N", null, true) }));
        neighbourBits.put(++bitpos, Arrays.asList(
                new Neighbour[] { new Neighbour("C", null, null), new Neighbour("C", IBond.Order.DOUBLE, false) }));
        smilesPatternBits.put("C:C:C:C:N:C", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + "1", ++bitpos);
        smilesPatternBits.put("O=C-C-C-C-C", ++bitpos);
        smilesPatternBits.put("C:C:C:C-C=O", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "1", ++bitpos);
        ringSetBits.put(hexRingMixNeighPrefix + "1", ++bitpos);
        smilesPatternBits.put("N-C-C-C:C", ++bitpos);
        smilesPatternBits.put("C-N-C-C-N-C", ++bitpos);
        elemCntBits.put("N5", ++bitpos);
        smilesPatternBits.put("O-C-C-C-C-C", ++bitpos);
        smilesPatternBits.put("O=C-C-N", ++bitpos);
        elemCntBits.put("C28", ++bitpos);
        smilesPatternBits.put("C:C:C-C-C-N", ++bitpos);
        smilesPatternBits.put("O-C-C-N-C", ++bitpos);
        smilesPatternBits.put("C:C:N:C:C:C", ++bitpos);
        smilesPatternBits.put("C-C:N:C", ++bitpos);
        smilesPatternBits.put("N:C:C-C", ++bitpos);
        smilesPatternBits.put("C-C-N-C-C:C", ++bitpos);
        smilesPatternBits.put("C-C-N-C:C:C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, true), new Neighbour("N", null, true) }));
        smilesPatternBits.put("C-N-C-C-C:C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixAny + "4", ++bitpos);
        smilesPatternBits.put("O=C-C-N-C", ++bitpos);
        smilesPatternBits.put("C-C=C-C", ++bitpos);
        ringSetBits.put(ringsetCountPerSet + "3", ++bitpos);
        ringSetBits.put(ringsetCountConnectedPrefix + "2", ++bitpos);
        smilesPatternBits.put("O-C:C:C:C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("N", null, false), new Neighbour("N", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("C", null, true), new Neighbour("C", null, true) }));
        smilesPatternBits.put("C-C-N-C-C-O", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("C", null, true), new Neighbour("N", null, true) }));
        smilesPatternBits.put("C-C-O-C-C", ++bitpos);
        elemCntBits.put("Cl1", ++bitpos);
        smilesPatternBits.put("C-C-N-C-C=O", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixCarbonOnly + "3", ++bitpos);
        smilesPatternBits.put("C:C-C-C-C-C", ++bitpos);
        smilesPatternBits.put("N:C:N:C", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixAny + "2", ++bitpos);
        smilesPatternBits.put("O=C-N-C-C-N", ++bitpos);
        ringBits.put(ringNonaromDoubleBond + "1_1", ++bitpos);
        smilesPatternBits.put("C-O-C-C-C-C", ++bitpos);
        smilesPatternBits.put("C-C-C=C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(
                new Neighbour[] { new Neighbour("C", null, null), new Neighbour("N", IBond.Order.DOUBLE, false) }));
        smilesPatternBits.put("O-C-C:C", ++bitpos);
        smilesPatternBits.put("O=C-N-C-C=O", ++bitpos);
        smilesPatternBits.put("O-C-C-C-N", ++bitpos);
        smilesPatternBits.put("C:C-C-C-C-N", ++bitpos);
        smilesPatternBits.put("O-C-C:C:C:C", ++bitpos);
        smilesPatternBits.put("C-C-O-C:C", ++bitpos);
        elemCntBits.put("F1", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, true) }));
        smilesPatternBits.put("C-C-O-C:C:C", ++bitpos);
        smilesPatternBits.put("C-C-C-O-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("Cl", null, false) }));
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "1", ++bitpos);
        smilesPatternBits.put("O=C-C=C", ++bitpos);
        smilesPatternBits.put("C-C:C:C-C", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "1", ++bitpos);
        smilesPatternBits.put("O=C-N-C:C:C", ++bitpos);
        smilesPatternBits.put("O=C-C-C-C:C", ++bitpos);
        elemCntBits.put("O7", ++bitpos);
        smilesPatternBits.put("C-C-C:C-C", ++bitpos);
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                        new Neighbour("C", IBond.Order.DOUBLE, false), new Neighbour("C", IBond.Order.SINGLE, false),
                        new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("O-C:C:C-C", ++bitpos);
        smilesPatternBits.put("C:C:C:C:C-Cl", ++bitpos);
        elemCntBits.put("C32", ++bitpos);
        smilesPatternBits.put("C-C-C=C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("N", null, true), new Neighbour("N", null, false) }));
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                        new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("C", IBond.Order.SINGLE, false),
                        new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("C-N-C-C-C-O", ++bitpos);
        smilesPatternBits.put("N-C:C-C", ++bitpos);
        elemCntBits.put("N6", ++bitpos);
        smilesPatternBits.put("C-C:C:N:C:C", ++bitpos);
        smilesPatternBits.put("O-C:C-C", ++bitpos);
        smilesPatternBits.put("C-C:C:C:C-C", ++bitpos);
        smilesPatternBits.put("N:C:C:N", ++bitpos);
        smilesPatternBits.put("O-C-C-N-C=O", ++bitpos);
        smilesPatternBits.put("O-C-C-O", ++bitpos);
        smilesPatternBits.put("C:C-C:C", ++bitpos);
        smilesPatternBits.put("C:N:C:C:C-C", ++bitpos);
        smilesPatternBits.put("C:C-C:C:C:C", ++bitpos);
        smilesPatternBits.put("N:C:C:C:C-C", ++bitpos);
        smilesPatternBits.put("C:C:C:N:C-C", ++bitpos);
        smilesPatternBits.put("N-C:N:C", ++bitpos);
        smilesPatternBits.put("N-C-N-C", ++bitpos);
        atomPairBits.put("N-N", ++bitpos);
        atomPairBits.put("S-O", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, true), new Neighbour("C", null, true) }));
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "1", ++bitpos);
        smilesPatternBits.put("N-C-C-O-C", ++bitpos);
        smilesPatternBits.put("O=C-C=C-C", ++bitpos);
        smilesPatternBits.put("N-C-C-C-N", ++bitpos);
        smilesPatternBits.put("O-C-C-O-C", ++bitpos);
        smilesPatternBits.put("C-C-C-N-C:C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("S", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, false) }));
        smilesPatternBits.put("C-N-C-C-C-N", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("S", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("O", IBond.Order.DOUBLE, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("S", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("S", null, false), new Neighbour("C", null, true) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                new Neighbour("N", IBond.Order.DOUBLE, false), new Neighbour("N", IBond.Order.SINGLE, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("S", null, false),
                new Neighbour("C", null, false), new Neighbour("O", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("S", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("C-O-C-C:C", ++bitpos);
        smilesPatternBits.put("N-C:C:C:C-C", ++bitpos);
        smilesPatternBits.put("C=C-C-C-C", ++bitpos);
        smilesPatternBits.put("C:C:C-C:C:C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + "2", ++bitpos);
        smilesPatternBits.put("C:C-C=C", ++bitpos);
        smilesPatternBits.put("C=C-C:C", ++bitpos);
        smilesPatternBits.put("N-C-N-C-C", ++bitpos);
        smilesPatternBits.put("C:C:C-C-O-C", ++bitpos);
        smilesPatternBits.put("O-C-C-C-C-N", ++bitpos);
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, false), new Neighbour("C", null, false),
                        new Neighbour("C", null, false), new Neighbour("C", null, false),
                        new Neighbour("C", null, false) }));
        smilesPatternBits.put("C-O-C:C:C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                new Neighbour("N", IBond.Order.DOUBLE, false), new Neighbour("C", IBond.Order.SINGLE, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, null),
                new Neighbour("C", IBond.Order.DOUBLE, false), new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("S-C:C:C:C:C", ++bitpos);
        smilesPatternBits.put("C-C:C-N-C", ++bitpos);
        smilesPatternBits.put("N-C-C=C", ++bitpos);
        smilesPatternBits.put("C-O-C-C-N-C", ++bitpos);
        smilesPatternBits.put("N:C:N:C:C:C", ++bitpos);
        smilesPatternBits.put("O=C-C-C-N", ++bitpos);
        ringSetBits.put(ringsetCountHexPrefix + "3", ++bitpos);
        smilesPatternBits.put("O-C-C-C-O", ++bitpos);
        smilesPatternBits.put("C:C:C:C-C=C", ++bitpos);
        smilesPatternBits.put("C:C:N:C:N:C", ++bitpos);
        smilesPatternBits.put("C-C-C:C:N:C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + "4", ++bitpos);
        smilesPatternBits.put("C-C:C-O-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("N", null, true), new Neighbour("N", null, false) }));
        smilesPatternBits.put("C:C-C-N-C=O", ++bitpos);
        smilesPatternBits.put("C:C:C-C=C-C", ++bitpos);
        smilesPatternBits.put("N:C-C:C", ++bitpos);
        smilesPatternBits.put("C-C-C-C=C-C", ++bitpos);
        smilesPatternBits.put("C:C-C-C-C-O", ++bitpos);
        atomPairBits.put("N-O", ++bitpos);
        smilesPatternBits.put("N:C:N-C", ++bitpos);
        smilesPatternBits.put("C-N-C-N-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                new Neighbour("C", IBond.Order.DOUBLE, false), new Neighbour("N", IBond.Order.SINGLE, false) }));
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "1", ++bitpos);
        smilesPatternBits.put("O=C-C-C-N-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("F", null, false) }));
        smilesPatternBits.put("N-C:C:C-C", ++bitpos);
        smilesPatternBits.put("N-C:C:C:N", ++bitpos);
        smilesPatternBits.put("F-C:C:C:C:C", ++bitpos);
        atomPairBits.put("S-N", ++bitpos);
        smilesPatternBits.put("C-C-O-C-C-O", ++bitpos);
        smilesPatternBits.put("C-C:C-C-C-C", ++bitpos);
        smilesPatternBits.put("O=C-C-C-C-N", ++bitpos);
        smilesPatternBits.put("C-C-C-C-C=C", ++bitpos);
        smilesPatternBits.put("N:C-C:C:C:C", ++bitpos);
        smilesPatternBits.put("C-O-C-C-C-O", ++bitpos);
        smilesPatternBits.put("N-C-O-C-C", ++bitpos);
        smilesPatternBits.put("O-C-C=C", ++bitpos);
        smilesPatternBits.put("C-C-S-C", ++bitpos);
        smilesPatternBits.put("C:C-C-C:C", ++bitpos);
        smilesPatternBits.put("C:C:C-C-C:C", ++bitpos);
        smilesPatternBits.put("C-C-C-O-C:C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("S", null, false),
                new Neighbour("N", null, false), new Neighbour("O", null, false) }));
        smilesPatternBits.put("C-C-C:C:C-O", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("S", null, false),
                new Neighbour("N", null, false), new Neighbour("O", null, false), new Neighbour("O", null, false) }));
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, false), new Neighbour("C", null, false),
                        new Neighbour("C", null, false), new Neighbour("C", null, false),
                        new Neighbour("O", null, false) }));
        smilesPatternBits.put("C-C:C:C-C-C", ++bitpos);
        smilesPatternBits.put("C-N-C-N-C-C", ++bitpos);
        smilesPatternBits.put("O-C-C-C=O", ++bitpos);
        smilesPatternBits.put("C:N:C:C:C:N", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + "2", ++bitpos);
        smilesPatternBits.put("O=S-C:C:C:C", ++bitpos);
        smilesPatternBits.put("C:N:C:C:C-N", ++bitpos);
        ringSetBits.put(ringsetCountPerSet + "4", ++bitpos);
        smilesPatternBits.put("C-O-C-C-C-N", ++bitpos);
        smilesPatternBits.put("C-C-C-N-C-N", ++bitpos);
        smilesPatternBits.put("N-C=N-C", ++bitpos);
        smilesPatternBits.put("O-C:C-O", ++bitpos);
        smilesPatternBits.put("N:C:C:C-C-C", ++bitpos);
        smilesPatternBits.put("C:N:C:C:N-C", ++bitpos);
        smilesPatternBits.put("C-N:C:N:C:C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "1", ++bitpos);
        smilesPatternBits.put("O-C-C=O", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, false),
                new Neighbour("C", null, false), new Neighbour("O", null, false) }));
        smilesPatternBits.put("N-C-C:C-C", ++bitpos);
        smilesPatternBits.put("O-C:C-C-C", ++bitpos);
        smilesPatternBits.put("C-C-C-O-C=O", ++bitpos);
        smilesPatternBits.put("N:C:N:C:C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("N", null, true), new Neighbour("N", null, true) }));
        smilesPatternBits.put("C-C:C-N-C-C", ++bitpos);
        smilesPatternBits.put("C:C-N-C-C-N", ++bitpos);
        smilesPatternBits.put("C-N-C-C=C-C", ++bitpos);
        smilesPatternBits.put("C-C-O-C-C-N", ++bitpos);
        smilesPatternBits.put("C:N:C:N:C-C", ++bitpos);
        smilesPatternBits.put("O-C-C-C-C=O", ++bitpos);
        smilesPatternBits.put("C-C-C-C:C:N", ++bitpos);
        smilesPatternBits.put("O-C-C=C-C", ++bitpos);
        smilesPatternBits.put("N=C-N-C", ++bitpos);
        smilesPatternBits.put("O=C-C:C-C", ++bitpos);
        smilesPatternBits.put("O=C-C-C:C", ++bitpos);
        smilesPatternBits.put("C:C:C-C:N:C", ++bitpos);
        smilesPatternBits.put("O-C:C-O-C", ++bitpos);
        smilesPatternBits.put("N-C:N:C:C:C", ++bitpos);
        smilesPatternBits.put("C:C:C:N-C-C", ++bitpos);
        smilesPatternBits.put("C-N-C-C:C-C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixAny + "5", ++bitpos);
        ringSetBits.put(ringsetCountTotalPrefix + "2", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("F", null, false) }));
        smilesPatternBits.put("C-C-O-C-N-C", ++bitpos);
        smilesPatternBits.put("O=C-C-C:C:C", ++bitpos);
        smilesPatternBits.put("N-C:N:C:N:C", ++bitpos);
        elemCntBits.put("F3", ++bitpos);
        elemCntBits.put("Cl2", ++bitpos);
        smilesPatternBits.put("N-C:C-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("F", null, false), new Neighbour("F", null, false) }));
        smilesPatternBits.put("C-N:C:C:C:C", ++bitpos);
        smilesPatternBits.put("C:C-N-C-C:C", ++bitpos);
        smilesPatternBits.put("N-S-C:C", ++bitpos);
        smilesPatternBits.put("O-C-C-C:C", ++bitpos);
        smilesPatternBits.put("C-C-N-C:N:C", ++bitpos);
        elemCntBits.put("S2", ++bitpos);
        smilesPatternBits.put("C:N:C:N-C-C", ++bitpos);
        smilesPatternBits.put("N-S-C:C:C:C", ++bitpos);
        ringSetBits.put(hexRingNonAromNeighPrefix + "1", ++bitpos);
        smilesPatternBits.put("C-N-C:C:C-C", ++bitpos);
        smilesPatternBits.put("O-C-O-C", ++bitpos);
        smilesPatternBits.put("C-C:C-O-C-C", ++bitpos);
        smilesPatternBits.put("C:C:C-C-C-O", ++bitpos);
        smilesPatternBits.put("C-C=C-N-C-C", ++bitpos);
        smilesPatternBits.put("C-C:C-C:C", ++bitpos);
        smilesPatternBits.put("C-S-C-C-C", ++bitpos);
        smilesPatternBits.put("C-C:C-C:C:C", ++bitpos);
        smilesPatternBits.put("N-C-N-C:C", ++bitpos);
        smilesPatternBits.put("N:C:C:C-N-C", ++bitpos);
        smilesPatternBits.put("C:C:N:C-N-C", ++bitpos);
        smilesPatternBits.put("C-C-N:C:C:N", ++bitpos);
        ringBits.put(ringNonaromDoubleBond + "2_1", ++bitpos);
        smilesPatternBits.put("C-C-O-C-C:C", ++bitpos);
        smilesPatternBits.put("C=C-C-N-C-C", ++bitpos);
        smilesPatternBits.put("N-C-C-C:C-C", ++bitpos);
        smilesPatternBits.put("O=C-N-C:C-C", ++bitpos);
        smilesPatternBits.put("C-C-C:C-O-C", ++bitpos);
        smilesPatternBits.put("N-C:C:N", ++bitpos);
        smilesPatternBits.put("O=C-C:C:C-C", ++bitpos);
        smilesPatternBits.put("Cl-C:C:C:C-C", ++bitpos);
        smilesPatternBits.put("C:C:N:C-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, false),
                new Neighbour("C", null, false), new Neighbour("N", null, false) }));
        smilesPatternBits.put("C-O-C-C=C", ++bitpos);
        smilesPatternBits.put("N:C:C:C:N-C", ++bitpos);
        smilesPatternBits.put("N-C-C-C:C:N", ++bitpos);
        smilesPatternBits.put("O=C-N-C=O", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "2", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "1", ++bitpos);
        smilesPatternBits.put("C-C-C-C:N:C", ++bitpos);
        smilesPatternBits.put("C:C:C-C-C=C", ++bitpos);
        smilesPatternBits.put("N=C-N-C-C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "2", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "1", ++bitpos);
        smilesPatternBits.put("O=C-C:C-N", ++bitpos);
        smilesPatternBits.put("O=C-C-O-C", ++bitpos);
        smilesPatternBits.put("O-C:C:C-O", ++bitpos);
        smilesPatternBits.put("C-S-C:C", ++bitpos);
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                        new Neighbour("C", IBond.Order.DOUBLE, false), new Neighbour("N", IBond.Order.SINGLE, false),
                        new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("O=C-C:N", ++bitpos);
        smilesPatternBits.put("S:C:C:C", ++bitpos);
        ringSetBits.put(hexRingMixNeighPrefix + "2", ++bitpos);
        elemCntBits.put("O12", ++bitpos);
        smilesPatternBits.put("C-O-C-O-C", ++bitpos);
        smilesPatternBits.put("O=C-C=C-N", ++bitpos);
        ringSetBits.put(ringsetCountPentPrefix + "2", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + "3", ++bitpos);
        smilesPatternBits.put("C:S:C-C", ++bitpos);
        smilesPatternBits.put("O=C-N-C-N", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + "2", ++bitpos);
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                        new Neighbour("N", IBond.Order.DOUBLE, false), new Neighbour("N", IBond.Order.SINGLE, false),
                        new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("Cl-C:C:C-C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixCarbonOnly + "4", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("N", null, false), new Neighbour("N", null, true), new Neighbour("N", null, true) }));
        carbonTrails.put(cTrailPrefix + "14", ++bitpos);
        elemCntBits.put("Br1", ++bitpos);
        smilesPatternBits.put("O-C-O-C-C", ++bitpos);
        smilesPatternBits.put("O=C-C:C-O", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("S", null, false), new Neighbour("N", null, false) }));
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "2", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(
                new Neighbour[] { new Neighbour("C", null, null), new Neighbour("N", IBond.Order.TRIPLE, false) }));
        smilesPatternBits.put("O-C-C-C=C", ++bitpos);
        ringSetBits.put(ringsetCountConnectedPrefix + "3", ++bitpos);
        smilesPatternBits.put("N=C-C-C", ++bitpos);
        smilesPatternBits.put("C:C:N:N:C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                new Neighbour("N", IBond.Order.TRIPLE, false), new Neighbour("C", IBond.Order.SINGLE, false) }));
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                        new Neighbour("N", IBond.Order.DOUBLE, false), new Neighbour("C", IBond.Order.SINGLE, false),
                        new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("O-C=C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("N", null, true), new Neighbour("N", null, true) }));
        smilesPatternBits.put("N:C:N:C:N:C", ++bitpos);
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, false), new Neighbour("C", null, false),
                        new Neighbour("C", null, false), new Neighbour("C", null, false),
                        new Neighbour("N", null, false) }));
        ringBits.put("5" + ringPrefixRing + ringPrefixAny + "3", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(
                new Neighbour[] { new Neighbour("N", null, null), new Neighbour("O", IBond.Order.DOUBLE, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("O", null, false),
                new Neighbour("C", null, true), new Neighbour("C", null, true) }));
        carbonTrails.put(cTrailPrefix + "15", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixAny + "1", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixNonArom + "1", ++bitpos);
        smilesPatternBits.put("C=C-C=C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, false),
                new Neighbour("O", null, false), new Neighbour("O", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("O", IBond.Order.SINGLE, false) }));
        ringSetBits.put(hexRingNonAromNeighPrefix + "2", ++bitpos);
        smilesPatternBits.put("C-C=C-C=C", ++bitpos);
        smilesPatternBits.put("O-C-C:C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("C", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("Cl-C:C-C", ++bitpos);
        smilesPatternBits.put("O=N-C:C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("Br", null, false), new Neighbour("C", null, true) }));
        ringSetBits.put(ringsetCountHexPrefix + "4", ++bitpos);
        ringSetBits.put(ringsetCountPerSet + "5", ++bitpos);
        carbonTrails.put(cTrailPrefix + "16", ++bitpos);
        smilesPatternBits.put("S:C:C:N", ++bitpos);
        ringSetBits.put(ringsetPairPrefix + "5_5", ++bitpos);
        smilesPatternBits.put("N-N-C-C", ++bitpos);
        smilesPatternBits.put("C:C-N-C:C", ++bitpos);
        smilesPatternBits.put("C-C-S-C-C", ++bitpos);
        elemCntBits.put("P1", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "2", ++bitpos);
        smilesPatternBits.put("O-C:C:C-N", ++bitpos);
        smilesPatternBits.put("O-C:C-N", ++bitpos);
        atomPairBits.put("P-O", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("P", null, false),
                new Neighbour("O", null, false), new Neighbour("O", null, false) }));
        smilesPatternBits.put("C-O-C=C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(
                new Neighbour[] { new Neighbour("P", null, null), new Neighbour("O", IBond.Order.DOUBLE, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("P", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("O", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("C:N-C:C", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "1", ++bitpos);
        smilesPatternBits.put("O=C-N-N", ++bitpos);
        ringBits.put(ringNonaromDoubleBond + "1_2", ++bitpos);
        smilesPatternBits.put("O=C-O-C-C-O", ++bitpos);
        smilesPatternBits.put("O=C-C-C=O", ++bitpos);
        smilesPatternBits.put("N=C-C=C", ++bitpos);
        smilesPatternBits.put("N-C-S-C", ++bitpos);
        ringSetBits.put(hexRingAromNeighPrefix + "1", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "3", ++bitpos);
        smilesPatternBits.put("O-N-C-C", ++bitpos);
        carbonTrails.put(cTrailPrefix + "17", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixAny + "1", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixNonArom + "1", ++bitpos);
        smilesPatternBits.put("S-C:C-C", ++bitpos);
        smilesPatternBits.put("O=S-C-C", ++bitpos);
        smilesPatternBits.put("C=N-N-C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + "5", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixCarbonOnly + "1", ++bitpos);
        smilesPatternBits.put("O-C:C:N", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "2", ++bitpos);
        smilesPatternBits.put("C-C=N-N-C", ++bitpos);
        ringLayout.put("6DLM", ++bitpos);
        smilesPatternBits.put("N-C:C-N", ++bitpos);
        smilesPatternBits.put("S-C:N:C", ++bitpos);
        elemCntBits.put("O16", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + "4", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(
                new Neighbour[] { new Neighbour("C", null, null), new Neighbour("C", IBond.Order.TRIPLE, false) }));
        smilesPatternBits.put("N-C:C:C-N", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                new Neighbour("C", IBond.Order.TRIPLE, false), new Neighbour("C", IBond.Order.SINGLE, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("O", null, false),
                new Neighbour("C", null, false), new Neighbour("P", null, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("O", null, false),
                new Neighbour("C", null, false), new Neighbour("N", null, false) }));
        smilesPatternBits.put("N-C:N:N", ++bitpos);
        smilesPatternBits.put("Cl-C:C-Cl", ++bitpos);
        smilesPatternBits.put("O=C-O-C:C", ++bitpos);
        smilesPatternBits.put("O-C:N:C:N:C", ++bitpos);
        carbonTrails.put(cTrailPrefix + "18", ++bitpos);
        ringSetBits.put(ringsetCountConnectedPrefix + "4", ++bitpos);
        smilesPatternBits.put("O-C-C:C-O", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "2", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("Cl", null, false) }));
        elemCntBits.put("I1", ++bitpos);
        ringLayout.put("5SLM", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "2", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + "3", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "1", ++bitpos);
        carbonTrails.put(cTrailPrefix + "19", ++bitpos);
        smilesPatternBits.put("N#C-C-C", ++bitpos);
        ringSetBits.put(ringsetCountPerSet + "6", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("S", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("O", IBond.Order.SINGLE, false) }));
        smilesPatternBits.put("N-N-C-N", ++bitpos);
        smilesPatternBits.put("Br-C:C:C-C", ++bitpos);
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("S", null, null),
                        new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("O", IBond.Order.DOUBLE, false),
                        new Neighbour("O", IBond.Order.SINGLE, false) }));
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("P", null, false) }));
        smilesPatternBits.put("S-C:C-N", ++bitpos);
        smilesPatternBits.put("Cl-C:C-O", ++bitpos);
        smilesPatternBits.put("S-C=N-C", ++bitpos);
        carbonTrails.put(cTrailPrefix + "20", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixAny + "4", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixCarbonOnly + "5", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(
                new Neighbour[] { new Neighbour("C", null, null), new Neighbour("S", IBond.Order.DOUBLE, false) }));
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "4", ++bitpos);
        ringSetBits.put(ringsetPairPrefix + "4_6", ++bitpos);
        smilesPatternBits.put("O=C-C-O-C=O", ++bitpos);
        ringLayout.put("5DLM", ++bitpos);
        smilesPatternBits.put("C-C-C#C", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "2", ++bitpos);
        elemCntBits.put("Na1", ++bitpos);
        smilesPatternBits.put("N-N-C:C", ++bitpos);
        atomPairBits.put("S-S", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + "3", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, true), new Neighbour("I", null, false) }));
        smilesPatternBits.put("Cl-C:C-O-C", ++bitpos);
        smilesPatternBits.put("S=C-N-C", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "3", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + "5", ++bitpos);
        ringSetBits.put(ringsetPairPrefix + "3_6", ++bitpos);
        ringLayout.put("6DLL", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(
                new Neighbour[] { new Neighbour("N", null, null), new Neighbour("N", IBond.Order.DOUBLE, false) }));
        ringSetBits.put(ringsetCountHexPrefix + "5", ++bitpos);
        ringLayout.put("6SLM", ++bitpos);
        elemCntBits.put("P2", ++bitpos);
        smilesPatternBits.put("Cl-C:C-C=O", ++bitpos);
        ringLayout.put("6SHH", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "1", ++bitpos);
        ringSetBits.put(ringsetCountPentPrefix + "3", ++bitpos);
        ringSetBits.put(hexRingAromNeighPrefix + "2", ++bitpos);
        smilesPatternBits.put("N#C-C-C-C", ++bitpos);
        elemCntBits.put("Br2", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "2", ++bitpos);
        ringSetBits.put(ringsetPairPrefix + "3_5", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "3", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "1", ++bitpos);
        smilesPatternBits.put("N#C-C=C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("Br", null, false), new Neighbour("C", null, false) }));
        smilesPatternBits.put("S-C:C:C-N", ++bitpos);
        smilesPatternBits.put("O-C-C=N", ++bitpos);
        smilesPatternBits.put("Cl-C-C-N-C", ++bitpos);
        smilesPatternBits.put("Cl-C-C-C", ++bitpos);
        smilesPatternBits.put("C-N-C-N-C-N", ++bitpos);
        ringSetBits.put(ringsetPairPrefix + "4_5", ++bitpos);
        smilesPatternBits.put("Br-C:C-C", ++bitpos);
        atomPairBits.put("O-O", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "3", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "3", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "1", ++bitpos);
        ringBits.put(ringNonaromDoubleBond + "2_2", ++bitpos);
        smilesPatternBits.put("O-S-C:C", ++bitpos);
        smilesPatternBits.put("O=N-C:C-N", ++bitpos);
        smilesPatternBits.put("C-O-C-O-C-O", ++bitpos);
        smilesPatternBits.put("C-N-C=N-C-N", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + "4", ++bitpos);
        smilesPatternBits.put("S-C:C-O", ++bitpos);
        ringSetBits.put(ringsetCountTotalPrefix + "3", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixAny + "2", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixNonArom + "2", ++bitpos);
        smilesPatternBits.put("N-C:O:C", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixAny + "5", ++bitpos);
        ringSetBits.put(ringsetCountPerSet + "8", ++bitpos);
        ringLayout.put("5SHH", ++bitpos);
        atomPairBits.put("P-N", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "3", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "1", ++bitpos);
        smilesPatternBits.put("Cl-C-C=O", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("P", null, false),
                new Neighbour("C", null, false), new Neighbour("C", null, false) }));
        smilesPatternBits.put("Cl-C-C-C-C", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "4", ++bitpos);
        smilesPatternBits.put("N-C-C:N:C-N", ++bitpos);
        elemCntBits.put("B1", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + "4", ++bitpos);
        smilesPatternBits.put("S-C-S-C", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "5", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "4", ++bitpos);
        elemCntBits.put("Si1", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("Br", null, false), new Neighbour("C", null, false), new Neighbour("C", null, false) }));
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "3", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "2", ++bitpos);
        smilesPatternBits.put("C-N-C-C:C:O", ++bitpos);
        elemCntBits.put("Pt1", ++bitpos);
        smilesPatternBits.put("O:C:C:C-C-N", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("Si", null, false) }));
        smilesPatternBits.put("C-N-C-C-C:O", ++bitpos);
        ringBits.put(ringNonaromDoubleBond + "1_3", ++bitpos);
        smilesPatternBits.put("O=N-C:C-O", ++bitpos);
        smilesPatternBits.put("Br-C-C-C", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("C", null, false), new Neighbour("I", null, false) }));
        atomPairBits.put("B-O", ++bitpos);
        smilesPatternBits.put("N-C-C:N:N-C", ++bitpos);
        ringSetBits.put(ringsetCountHexPrefix + "6", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("N", null, true), new Neighbour("N", null, true), new Neighbour("Cl", null, false) }));
        smilesPatternBits.put("O=S-C-N", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "2", ++bitpos);
        ringLayout.put("5DHM", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("Br", null, false), new Neighbour("N", null, true) }));
        smilesPatternBits.put("Br-C-C=O", ++bitpos);
        smilesPatternBits.put("Cl-C-C-O", ++bitpos);
        atomPairBits.put("Si-O", ++bitpos);
        ringSetBits.put(ringsetCountPentPrefix + "4", ++bitpos);
        smilesPatternBits.put("C:C:N:C-C:O", ++bitpos);
        smilesPatternBits.put("O=C-N-C-C:O", ++bitpos);
        smilesPatternBits.put("C-C-C-N:C-O", ++bitpos);
        smilesPatternBits.put("C:N:C-C-C:N", ++bitpos);
        smilesPatternBits.put("C=N-C=N-C-N", ++bitpos);
        smilesPatternBits.put("C-N-C-O-C=O", ++bitpos);
        smilesPatternBits.put("C:N:C-C-N=C", ++bitpos);
        ringLayout.put("6DHM", ++bitpos);
        smilesPatternBits.put("O-C=C-N-C-C", ++bitpos);
        smilesPatternBits.put("C:C-N-C-C:O", ++bitpos);
        ringLayout.put("5DLL", ++bitpos);
        smilesPatternBits.put("C:N:C:C-C:O", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, false),
                new Neighbour("Cl", null, false), new Neighbour("Cl", null, false) }));
        ringSetBits.put(ringsetCountPerSet + "10", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + "5", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "5", ++bitpos);
        smilesPatternBits.put("C:C:C:C=C-C", ++bitpos);
        ringLayout.put("5DHH", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "4", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + "5", ++bitpos);
        ringLayout.put("6DHH", ++bitpos);
        smilesPatternBits.put("O=C-C:N:C:O", ++bitpos);
        smilesPatternBits.put("O:C:C:C-N-C", ++bitpos);
        smilesPatternBits.put("C-C-N:C:N-C", ++bitpos);
        smilesPatternBits.put("N:C:N:N-C:C", ++bitpos);
        smilesPatternBits.put("O-C-C-C-C:O", ++bitpos);
        smilesPatternBits.put("O=C-N-C=C-O", ++bitpos);
        smilesPatternBits.put("N:C:N:C-N:C", ++bitpos);
        smilesPatternBits.put("C:C:N-C-C:N", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "4", ++bitpos);
        ringLayout.put("5SHM", ++bitpos);
        smilesPatternBits.put("O:C-C-C-C-N", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "3", ++bitpos);
        smilesPatternBits.put("S-S-C:C", ++bitpos);
        smilesPatternBits.put("C-C-O-C-C:O", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixCarbonOnly + "2", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "2", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "3", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "5", ++bitpos);
        smilesPatternBits.put("C:N:C:C:C:O", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "4", ++bitpos);
        smilesPatternBits.put("O:C:C:C:C:N", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixAny + "2", ++bitpos);
        elemCntBits.put("K1", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixNonArom + "2", ++bitpos);
        smilesPatternBits.put("N:C:N:N:N:C", ++bitpos);
        smilesPatternBits.put("O=C-C-C-C:O", ++bitpos);
        ringSetBits.put(ringsetCountTotalPrefix + "4", ++bitpos);
        elemCntBits.put("Se1", ++bitpos);
        smilesPatternBits.put("N-C-C-C:C:O", ++bitpos);
        smilesPatternBits.put("O:C-N-C:C:C", ++bitpos);
        ringLayout.put("5SLL", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "4", ++bitpos);
        atomPairBits.put("Cl-O", ++bitpos);
        ringBits.put("3" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "2", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "5", ++bitpos);
        ringSetBits.put(ringsetPairPrefix + "3_3", ++bitpos);
        atomPairBits.put("F-S", ++bitpos);
        smilesPatternBits.put("C:S:C:C-C:O", ++bitpos);
        smilesPatternBits.put("N:N:C-C-C:N", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "4", ++bitpos);
        atomPairBits.put("P-S", ++bitpos);
        elemCntBits.put("Li1", ++bitpos);
        ringSetBits.put(ringsetCountPerSet + "30", ++bitpos);
        elemCntBits.put("Tc1", ++bitpos);
        elemCntBits.put("Au1", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, false),
                new Neighbour("O", null, false), new Neighbour("O", null, true) }));
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "5", ++bitpos);
        ringLayout.put("5DLH", ++bitpos);
        ringLayout.put("6DLH", ++bitpos);
        elemCntBits.put("Fe1", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixNitro + "5", ++bitpos);
        ringLayout.put("5SLH", ++bitpos);
        atomPairBits.put("B-N", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixNonArom + ringPrefixHetero + "2", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "3", ++bitpos);
        smilesPatternBits.put("Br-C-C-C:C", ++bitpos);
        ringSetBits.put(ringsetPairPrefix + "3_4", ++bitpos);
        elemCntBits.put("Cu1", ++bitpos);
        ringLayout.put("6SLL", ++bitpos);
        ringLayout.put("6SHM", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("N", null, false),
                new Neighbour("C", null, true), new Neighbour("C", null, true), new Neighbour("C", null, true) }));
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "5", ++bitpos);
        atomPairBits.put("B-F", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "2", ++bitpos);
        smilesPatternBits.put("Cl-C-C-Cl", ++bitpos);
        ringBits.put(ringNonaromDoubleBond + "2_3", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixCarbonOnly + "3", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "4", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "5", ++bitpos);
        elemCntBits.put("Ru1", ++bitpos);
        elemCntBits.put(otherElem, ++bitpos);
        elemCntBits.put("Mn1", ++bitpos);
        elemCntBits.put("Zn1", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "4", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixCarbonOnly + "4", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixArom + ringPrefixCarbonOnly + "5", ++bitpos);
        atomPairBits.put("P-F", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(
                new Neighbour[] { new Neighbour("N", null, null), new Neighbour("P", IBond.Order.DOUBLE, false) }));
        elemCntBits.put("Re1", ++bitpos);
        elemCntBits.put("Sn1", ++bitpos);
        elemCntBits.put("As1", ++bitpos);
        elemCntBits.put("Te1", ++bitpos);
        atomPairBits.put("P-Cl", ++bitpos);
        elemCntBits.put("Co1", ++bitpos);
        ringBits.put("4" + ringPrefixRing + ringPrefixNonArom + ringPrefixNitro + "2", ++bitpos);
        ringLayout.put("5SHL", ++bitpos);
        ringLayout.put("5DHL", ++bitpos);
        elemCntBits.put("Bi1", ++bitpos);
        elemCntBits.put("Ni1", ++bitpos);
        elemCntBits.put("Pd1", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNonArom + ringPrefixCarbonOnly + "5", ++bitpos);
        neighbourBits.put(++bitpos, Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("Cl", IBond.Order.SINGLE, false) }));
        elemCntBits.put("Sb1", ++bitpos);
        neighbourBits.put(++bitpos,
                Arrays.asList(new Neighbour[] { new Neighbour("C", null, null),
                        new Neighbour("O", IBond.Order.DOUBLE, false), new Neighbour("Cl", IBond.Order.SINGLE, false),
                        new Neighbour("C", IBond.Order.SINGLE, false) }));
        elemCntBits.put("Ca1", ++bitpos);
        elemCntBits.put("V1", ++bitpos);
        ringLayout.put("6SLH", ++bitpos);
        ringLayout.put("6DHL", ++bitpos);
        elemCntBits.put("Cr1", ++bitpos);
        elemCntBits.put("Hg1", ++bitpos);
        atomPairBits.put("Cl-N", ++bitpos);
        ringLayout.put("6SHL", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixArom + ringPrefixHetero + "5", ++bitpos);

        // new for SMARTS
        ringBits.put("5" + ringPrefixRing + ringPrefixHetero + "1", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixHetero + "2", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixHetero + "4", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixHetero + "1", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixHetero + "2", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixHetero + "4", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNitro + "1", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNitro + "2", ++bitpos);
        ringBits.put("5" + ringPrefixRing + ringPrefixNitro + "4", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNitro + "1", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNitro + "2", ++bitpos);
        ringBits.put("6" + ringPrefixRing + ringPrefixNitro + "4", ++bitpos);

        ringSizes.put(ringSize + "3", ++bitpos);
        ringSizes.put(ringSize + "4", ++bitpos);

        // groups of related and misleadingly rare smiles pattern
        bitpos += 1;
        smilesPatternBits.put("O:C:C-C:C:O", bitpos);
        smilesPatternBits.put("O:C-O-C-O-C", bitpos);
        smilesPatternBits.put("O:C-O-C-C=O", bitpos);
        smilesPatternBits.put("O:C-O-C-C-O", bitpos);
        smilesPatternBits.put("O:C-O-C-C:C", bitpos);
        smilesPatternBits.put("O:C-C-O-C-O", bitpos);
        smilesPatternBits.put("C:C-O-C:C:O", bitpos);
        bitpos += 1;
        smilesPatternBits.put("C:C=C-C-O-C", bitpos);
        smilesPatternBits.put("O:C:C:C:C=C", bitpos);
        smilesPatternBits.put("O=C-C:C=C-C", bitpos);
        smilesPatternBits.put("C-O-C:C=C-O", bitpos);
        smilesPatternBits.put("C:C:C=C-C-O", bitpos);
        smilesPatternBits.put("C=C:C-C-C-O", bitpos);
        smilesPatternBits.put("C-C:C:C=C-O", bitpos);
        smilesPatternBits.put("O:C-C=C-C:O", bitpos);
        smilesPatternBits.put("O=C-C=C:C-C", bitpos);
        smilesPatternBits.put("C-C=C:C-C-O", bitpos);
        smilesPatternBits.put("O-C:C:C=C-O", bitpos);
        bitpos += 1;
        smilesPatternBits.put("O-C-C-C=C=C", bitpos);
        smilesPatternBits.put("C=C=C-C-C=O", bitpos);
        smilesPatternBits.put("O-O-O-C-C-C", bitpos);
        smilesPatternBits.put("O-C-C=C=C-C", bitpos);
        bitpos += 1;
        smilesPatternBits.put("C-C=C=C=C-C", bitpos);
        smilesPatternBits.put("C=C=C-C-C=C", bitpos);
        smilesPatternBits.put("C-C=C-C=C=C", bitpos);
        smilesPatternBits.put("C-C=C=C-C=C", bitpos);
        bitpos += 1;
        smilesPatternBits.put("C-C=C-C:C=C", bitpos);
        smilesPatternBits.put("C:C-C=C=C=C", bitpos);
        smilesPatternBits.put("C:C:C-C:C=C", bitpos);
        smilesPatternBits.put("C:C-C=C:C=C", bitpos);
        smilesPatternBits.put("C:C-C-C:C=C", bitpos);
        smilesPatternBits.put("C:C=C-C-C:C", bitpos);
        smilesPatternBits.put("C-C=C=C-C:C", bitpos);
        smilesPatternBits.put("C=C-C=C:C-C", bitpos);
        smilesPatternBits.put("C=C:C:C:C-C", bitpos);
        smilesPatternBits.put("C-C=C:C=C-C", bitpos);
        smilesPatternBits.put("C:C-C=C:C-C", bitpos);
        smilesPatternBits.put("C-C=C:C-C=C", bitpos);
        smilesPatternBits.put("C-C-C=C:C=C", bitpos);
        smilesPatternBits.put("C=C:C:C-C=C", bitpos);
        smilesPatternBits.put("C-C:C-C=C:C", bitpos);
        smilesPatternBits.put("C-C=C:C-C-C", bitpos);
        smilesPatternBits.put("C-C:C=C-C-C", bitpos);
        smilesPatternBits.put("C-C:C:C=C-C", bitpos);
        smilesPatternBits.put("C-C-C-C=C:C", bitpos);
        bitpos += 1;
        smilesPatternBits.put("N-C-C:C=C-C", bitpos);
        smilesPatternBits.put("N=C:C-N-C-N", bitpos);
        smilesPatternBits.put("C=C:N:C-C-C", bitpos);
        smilesPatternBits.put("C=C:C-C-N-C", bitpos);
        smilesPatternBits.put("C-C-C:C=C-N", bitpos);
        smilesPatternBits.put("N-C-N=C:C-N", bitpos);
        smilesPatternBits.put("C=C:C:N:C-C", bitpos);
        smilesPatternBits.put("N-C-C=C:C-C", bitpos);
        smilesPatternBits.put("N-C-C-C=C:C", bitpos);
        smilesPatternBits.put("C-N-C:C=C-C", bitpos);
        smilesPatternBits.put("C=C-C:N:C=C", bitpos);
        smilesPatternBits.put("N-C-C=C:C:N", bitpos);
        smilesPatternBits.put("N-C:C=N-C-C", bitpos);
        bitpos += 1;
        smilesPatternBits.put("C-N:C-C:N-N", bitpos);
        smilesPatternBits.put("N-N-N:C:C:C", bitpos);
        smilesPatternBits.put("N-N-N:C-C-C", bitpos);
        smilesPatternBits.put("N-N:C-C-C-N", bitpos);
        smilesPatternBits.put("N:C-N:C-C:N", bitpos);
        smilesPatternBits.put("N:C-N-C:N-C", bitpos);
        smilesPatternBits.put("N:C-C:N-C-N", bitpos);
        smilesPatternBits.put("N:C-C-C:N-N", bitpos);
        bitpos += 1;
        smilesPatternBits.put("C:C-N:N-C:C", bitpos);
        smilesPatternBits.put("N:C:N:N:N-N", bitpos);
        smilesPatternBits.put("N-N:N:C:N:N", bitpos);
        smilesPatternBits.put("C-N:N:C:N-N", bitpos);
        smilesPatternBits.put("N:C:N:N-N-C", bitpos);
        smilesPatternBits.put("N:N-C-C-N:N", bitpos);
        smilesPatternBits.put("C-N:N:N-C-C", bitpos);
        smilesPatternBits.put("N:N:C-C:N-N", bitpos);
        smilesPatternBits.put("C-N:N-C-C-N", bitpos);
        smilesPatternBits.put("C:N:N:N:N:C", bitpos);
        smilesPatternBits.put("N-N:N:C-C-C", bitpos);
        smilesPatternBits.put("C:C:N:N:N:N", bitpos);
        smilesPatternBits.put("N-N:N:N:C-C", bitpos);
        smilesPatternBits.put("N:N-N-C-C-C", bitpos);
        smilesPatternBits.put("C-N-N:N:N:C", bitpos);
        smilesPatternBits.put("C:C:N:N:N-N", bitpos);
        smilesPatternBits.put("N-N:N-C-C-C", bitpos);
        smilesPatternBits.put("N:N:N:N:C-O", bitpos);
        bitpos += 1;
        smilesPatternBits.put("C-C-C-N=C=N", bitpos);
        smilesPatternBits.put("C-C-N=C=N-C", bitpos);
        smilesPatternBits.put("C=C=C-N-C-C", bitpos);
        smilesPatternBits.put("N=C-C=N-N-N", bitpos);
        bitpos += 1;
        smilesPatternBits.put("C-C:C:S:C=N", bitpos);
        smilesPatternBits.put("S:C:C:S:C:N", bitpos);
        smilesPatternBits.put("N:S:C:C-S-N", bitpos);
        smilesPatternBits.put("N:N:N:N-C-N", bitpos);
        smilesPatternBits.put("N-S-N-S-C:C", bitpos);
        smilesPatternBits.put("C:S:C:S:C-N", bitpos);
        smilesPatternBits.put("C-N-C:S:C:S", bitpos);
        smilesPatternBits.put("C:S:C:S:C:N", bitpos);
        smilesPatternBits.put("C-N=C:S:C:C", bitpos);
        bitpos += 1;
        smilesPatternBits.put("O-C=C-N-C:O", bitpos);
        smilesPatternBits.put("C-C-N-C=C:O", bitpos);
        smilesPatternBits.put("O=C-N-C=C:O", bitpos);
        smilesPatternBits.put("N-C:C:C=C-O", bitpos);
        smilesPatternBits.put("O-C-N-C=C:O", bitpos);
        smilesPatternBits.put("C:C:N:C=C-O", bitpos);
        smilesPatternBits.put("O:C-C=C:C:N", bitpos);
        smilesPatternBits.put("C=C-C-N-C:O", bitpos);
        bitpos += 1;
        smilesPatternBits.put("N-C:O:N:C-O", bitpos);
        smilesPatternBits.put("C:C-C:O:N-C", bitpos);
        smilesPatternBits.put("O:C:C:N-C-O", bitpos);
        smilesPatternBits.put("O:C:N-C-C-C", bitpos);
        smilesPatternBits.put("O-C:N:N:N-C", bitpos);
        smilesPatternBits.put("O:C-C-O-C:N", bitpos);
        smilesPatternBits.put("O:C-O-C-C-N", bitpos);
        smilesPatternBits.put("O=N-C:C:O:N", bitpos);
        smilesPatternBits.put("O-C-N-C:O:N", bitpos);
        smilesPatternBits.put("O:C-N-C-C:O", bitpos);
        smilesPatternBits.put("C-C-N:C:C:O", bitpos);
        bitpos += 1;
        smilesPatternBits.put("C=N-C-N=C-O", bitpos);
        smilesPatternBits.put("O=C-C-N=C=O", bitpos);
        smilesPatternBits.put("O-C-C-N=C=O", bitpos);
        bitpos += 1;
        smilesPatternBits.put("N:C:S:C:C:O", bitpos);
        smilesPatternBits.put("O:C-N-S-C:N", bitpos);
        smilesPatternBits.put("O:C:N:C:S:C", bitpos);
        bitpos += 1;
        smilesPatternBits.put("S:C-C:C-C:O", bitpos);
        smilesPatternBits.put("S-C:S:C:C:O", bitpos);
        smilesPatternBits.put("O:C-C-C:C:S", bitpos);
        smilesPatternBits.put("S:C-C-C-C:O", bitpos);
        smilesPatternBits.put("O-C:S:C-S=O", bitpos);
        bitpos += 1;
        smilesPatternBits.put("C-O-S=N-C-C", bitpos);
        smilesPatternBits.put("C-C-C-N=S-O", bitpos);

        /* Make a list (for reference) of elements for which counts are done is this fingerprint */
        Iterator<String> iterator = elemCntBits.keySet().iterator();

        while(iterator.hasNext())
        {
            String elemSymbol = iterator.next();
            elemSymbol = elemSymbol.replaceAll("[0-9]", "");
            elemFingerprinted.put(elemSymbol, elemSymbol);
        }
    }
}
