/*
 * Copyright (C) 2015-2017 Jakub Galgonek   galgonek@uochb.cas.cz
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
package cz.iocb.sachem.tautomers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;



public class InChI
{
    public static class Fragment
    {
        public String formula;
        public String connections;
        public String hydrogens;
        public String formalCharge;
        public String doubleBondStereo;
        public String tetrahedralStereo;
        String aux;


        void setPart(char name, String part)
        {
            switch(name)
            {
                case 'c':
                    connections = part;
                    break;
                case 'h':
                    hydrogens = part;
                    break;
                case 'q':
                    formalCharge = part;
                    break;
                case 'b':
                    doubleBondStereo = part;
                    break;
                case 't':
                    tetrahedralStereo = part;
                    break;
            }
        }
    }


    static private final ProcessBuilder pb;

    private String key = null;
    private String aux = null;
    private String value = null;


    static
    {
        String path = System.getProperty("inchi.path");
        pb = new ProcessBuilder(path, "-STDIO", "-Key", "-NoLabels", "-W0", "-SUU", "-RecMet");
    }


    public InChI(String molfile)
    {
        try
        {
            Process p = pb.start();
            BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
            BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            BufferedWriter output = new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));

            output.write(molfile);
            output.close();

            value = input.readLine();
            aux = input.readLine();
            key = input.readLine();

            // workaround
            if(value != null && value.contains("/r"))
            {
                value = "InChI=1/" + value.substring(value.indexOf("/r") + 2);
                aux = "AuxInfo=1/" + aux.substring(aux.indexOf("/R:") + 4);
            }

            input.close();
            error.close();
        }
        catch(IOException e)
        {
            e.printStackTrace();
        }
    }


    public List<Fragment> decompose()
    {
        LinkedList<Fragment> fragments = new LinkedList<Fragment>();

        if(!aux.contains("/N:"))
            return fragments;

        String[] components = value.split("/");

        String formulaComponent = components[1];
        String[] formulas = formulaComponent.split("\\.");

        for(String formula : formulas)
        {
            Pattern formulaPattern = Pattern.compile("^[0-9]+");
            Matcher match = formulaPattern.matcher(formula);

            if(!match.find())
            {
                Fragment fragment = new Fragment();
                fragment.formula = formula;
                fragments.add(fragment);
            }
            else
            {
                String countString = match.group();
                String value = formula.substring(countString.length());
                int count = Integer.valueOf(countString);

                for(int i = 0; i < count; i++)
                {
                    Fragment fragment = new Fragment();
                    fragment.formula = value;
                    fragments.add(fragment);
                }
            }
        }


        boolean hasIsomericPart = false;

        for(int i = 2; i < components.length; i++)
        {
            char name = components[i].charAt(0);

            if(name == 'p' || name == 'h' && hasIsomericPart)
                continue;

            if(name == 'i')
                hasIsomericPart = true;

            String[] parts = components[i].substring(1).split(";");
            int index = 0;

            for(String part : parts)
            {
                Pattern formulaPattern = Pattern.compile("^[0-9]+\\*");
                Matcher match = formulaPattern.matcher(part);

                if(!match.find())
                {
                    if(!part.isEmpty())
                        fragments.get(index).setPart(name, part);

                    index++;
                }
                else
                {
                    String countString = match.group();
                    String value = part.substring(countString.length());
                    int count = Integer.valueOf(countString.substring(0, countString.length() - 1));

                    for(int j = 0; j < count; j++)
                    {
                        if(!value.isEmpty())
                            fragments.get(index).setPart(name, value);

                        index++;
                    }
                }
            }
        }


        String[] auxComponents = aux.split("/");

        if(auxComponents[2].startsWith("N:"))
        {
            String[] values = auxComponents[2].substring(2).split(";");

            int index = 0;
            for(String value : values)
                fragments.get(index++).aux = value;
        }

        return fragments;
    }


    public String getKey()
    {
        return key;
    }


    public String getAuxInfo()
    {
        return aux;
    }


    public String getValue()
    {
        return value;
    }
}
