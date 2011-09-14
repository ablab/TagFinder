package ru.spbau.bioinf.tagfinder;

import java.util.HashMap;

public class Acids {

    public static HashMap<Character, Double> acids = new HashMap<Character, Double>();

    static {
        acids.put('G', 57.02146);
        acids.put('A', 71.03711);
        acids.put('S', 87.03203);
        acids.put('P', 97.05276);
        acids.put('V', 99.06841);
        acids.put('T', 101.04768);
        acids.put('C', 103.00919);
        acids.put('I', 113.08406);
        //acids.put('L', 113.08406);
        acids.put('N', 114.04293);
        acids.put('D', 115.02694);
        acids.put('R', 156.10111);
        acids.put('Q', 128.05858);
        acids.put('K', 128.09496);
        acids.put('E', 129.04259);
        acids.put('M', 131.04049);
        acids.put('H', 137.05891);
        acids.put('F', 147.06841);
        acids.put('U', 150.953636);
        acids.put('Y', 163.06333);
        acids.put('W', 186.07931);
        acids.put('O', 237.147727);
    }
}