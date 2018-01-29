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
package cz.iocb.sachem.shared;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import cz.iocb.sachem.search.SachemMoleculeBuilder;
import cz.iocb.sachem.search.SachemMoleculeBuilder.BondType;



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


    public MoleculeCounts(IAtomContainer iac, boolean isTarget) throws CDKException
    {
        for(IAtom atom : iac.atoms())
        {
            if(!(atom instanceof IPseudoAtom))
            {
                if(!atom.getSymbol().equals("H"))
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
            else if(isTarget)
            {
                molAtomCount++;
                molSCount++;
                molNCount++;
                molOCount++;
                molFCount++;
                molClCount++;
                molBrCount++;
                molICount++;
                molCCount++;
                molPCount++;
            }
        }


        bondLoop:
        for(IBond bond : iac.bonds())
        {
            for(IAtom atom : bond.atoms())
                if(atom.getSymbol().equals("H"))
                    continue bondLoop;


            BondType type = SachemMoleculeBuilder.getBondType(bond);

            if(type == BondType.AROMATIC)
                molAromaticBondCount++;
            else if(type == BondType.SINGLE)
                molSingleBondCount++;
            else if(type == BondType.DOUBLE)
                molDoubleBondCount++;
            else if(type == BondType.TRIPLE)
                molTripleBondCount++;

            if(isTarget)
            {
                if(type == BondType.SINGLE_OR_DOUBLE)
                {
                    molSingleBondCount++;
                    molDoubleBondCount++;
                }
                else if(type == BondType.SINGLE_OR_AROMATIC)
                {
                    molSingleBondCount++;
                    molAromaticBondCount++;
                }
                else if(type == BondType.DOUBLE_OR_AROMATIC)
                {
                    molDoubleBondCount++;
                    molAromaticBondCount++;
                }
                else if(type == BondType.ANY)
                {
                    molSingleBondCount++;
                    molDoubleBondCount++;
                    molTripleBondCount++;
                    molAromaticBondCount++;
                }
            }
        }
    }
}
