/*
 * Copyright (C) 2015-2017 Jakub Galgonek   galgonek@uochb.cas.cz
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
package cz.iocb.orchem.shared;

import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;



/**
 * Class to counting (particular) elements and bond orders in a given atom container
 */
public class MoleculeCounts
{
    public short molSingleBondCount;
    public short molDoubleBondCount;
    public short molTripleBondCount;
    public short molAromaticBondCount;
    public short molSCount;
    public short molOCount;
    public short molNCount;
    public short molFCount;
    public short molClCount;
    public short molBrCount;
    public short molICount;
    public short molCCount;
    public short molPCount;
    public short molAtomCount;


    public MoleculeCounts(IAtomContainer iac)
    {
        for(IAtom atom : iac.atoms())
        {
            molAtomCount++;

            if(atom.getSymbol().equals("S"))
                molSCount++;
            else if(atom.getSymbol().equals("N"))
                molNCount++;
            else if(atom.getSymbol().equals("O"))
                molOCount++;
            else if(atom.getSymbol().equals("F"))
                molFCount++;
            else if(atom.getSymbol().equals("Cl"))
                molClCount++;
            else if(atom.getSymbol().equals("Br"))
                molBrCount++;
            else if(atom.getSymbol().equals("I"))
                molICount++;
            else if(atom.getSymbol().equals("C"))
                molCCount++;
            else if(atom.getSymbol().equals("P"))
                molPCount++;
        }

        for(IBond bond : iac.bonds())
        {
            if(bond.getFlag(CDKConstants.ISAROMATIC))
                molAromaticBondCount++;
            else if(bond.getOrder() == Bond.Order.SINGLE)
                molSingleBondCount++;
            else if(bond.getOrder() == Bond.Order.DOUBLE)
                molDoubleBondCount++;
            else if(bond.getOrder() == Bond.Order.TRIPLE)
                molTripleBondCount++;
        }
    }
}
