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
package cz.iocb.orchem.tautomers;

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
        public String charge;
        public String fixedHydrogens;
        public String formalCharge;
        public int protonation;
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
                case 'f':
                    fixedHydrogens = part;
                    break;
                case 'q':
                    formalCharge = part;
                    break;
            }
        }
    }


    static private final ProcessBuilder pb;

    private String key = null;
    private String aux = null;
    private String value = null;
    private String fixedH = null;


    static
    {
        String path = System.getProperty("inchi.path");
        pb = new ProcessBuilder(path, "-STDIO", "-Key", "-NoLabels", "-W0", "-FixedH", "-RecMet");
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


            if(value != null)
            {
                int fixedHIdx = value.indexOf("/f");

                if(fixedHIdx != -1)
                {
                    fixedH = value.substring(fixedHIdx + 2);
                    value = value.substring(0, fixedHIdx);
                }
                else
                {
                    fixedH = "";
                }
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

        for(int i = 2; i < components.length; i++)
        {
            char name = components[i].charAt(0);

            if(name == 'p')
                continue;

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


        if(!fixedH.startsWith("/"))
        {
            String orgFormulaComponent = fixedH.substring(0, fixedH.indexOf('/', 0));
            String[] orgFormulas = orgFormulaComponent.split("\\.");

            int index = 0;

            for(String orgFormula : orgFormulas)
            {
                Pattern formulaPattern = Pattern.compile("^[0-9]+");
                Matcher match = formulaPattern.matcher(orgFormula);

                int count = 1;

                if(match.find())
                {
                    count = Integer.valueOf(match.group());
                    orgFormula = orgFormula.substring(match.group().length());
                }


                for(int i = 0; i < count; i++)
                {
                    if(index < fragments.size())
                    {
                        Fragment fragment = fragments.get(index);
                        fragment.protonation = countOfHydrogens(orgFormula) - countOfHydrogens(fragment.formula);
                        fragment.formula = orgFormula;
                    }
                    else
                    {
                        Fragment fragment = new Fragment();
                        fragment.formula = orgFormula;
                        fragments.add(fragment);
                    }

                    index++;
                }
            }
        }


        String[] fixedHcomponents = fixedH.split("/");

        for(int i = 1; i < fixedHcomponents.length; i++)
        {
            char name = fixedHcomponents[i].charAt(0);

            if(name == 'h')
                name = 'f';

            String[] parts = fixedHcomponents[i].substring(1).split(";");
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


    private int countOfHydrogens(String formula)
    {
        if(!formula.contains("H"))
            return 0;

        Pattern hPattern = Pattern.compile("H[0-9]+");
        Matcher match = hPattern.matcher(formula);

        if(match.find())
            return Integer.valueOf(match.group().substring(1));

        return 1;
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


    public String getFixedH()
    {
        return fixedH;
    }
}
