/*
 * Copyright (C) 2015-2017 Jakub Galgonek   galgonek@uochb.cas.cz
 * Copyright (C) 2011-2011 Mark Rijnbeek    markr@ebi.ac.uk
 * Copyright (C) 2008-2009 Federico Paoli
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
package cz.iocb.sachem.convert;

import java.io.IOException;
import java.io.StringWriter;
import java.util.Properties;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;



public class ConvertMolecule
{
    public static String smilesToMolfile(String smiles, boolean generateCoords, boolean useBondType4)
            throws CDKException, IOException
    {
        String cmolfile = null;

        if(smiles != null && !smiles.trim().equals(""))
        {
            SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
            sp.kekulise(!useBondType4);

            IAtomContainer molecule = sp.parseSmiles(smiles);

            // why are valencies being set?
            for(IAtom atom : molecule.atoms())
                atom.setValency(null);

            if(generateCoords)
            {
                StructureDiagramGenerator sdg = new StructureDiagramGenerator();
                sdg.setMolecule(molecule);
                sdg.generateCoordinates();
                molecule = sdg.getMolecule();
            }

            StringWriter out = new StringWriter();
            MDLV2000Writer mdlWriter = new MDLV2000Writer(out);

            if(useBondType4)
            {
                Properties prop = new Properties();
                prop.setProperty("WriteAromaticBondTypes", "true");
                PropertiesListener listener = new PropertiesListener(prop);
                mdlWriter.addChemObjectIOListener(listener);
                mdlWriter.customizeJob();
            }

            mdlWriter.setWriter(out);
            mdlWriter.write(molecule);
            mdlWriter.close();
            cmolfile = out.toString();
            out.close();
        }

        return cmolfile;
    }
}
