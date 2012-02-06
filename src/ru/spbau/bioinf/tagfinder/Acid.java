package ru.spbau.bioinf.tagfinder;

import java.util.HashMap;

public enum Acid {

    G(57.02146, 2, 3, 1, 1, 0),
    A(71.03711, 3, 5, 1, 1, 0),
    S(87.03203, 3, 5, 1, 2, 0),
    P(97.05276, 5, 7, 1, 1, 0),
    V(99.06841, 5, 9, 1, 1, 0),
    T(101.04768, 4, 7, 1, 2, 0),
    C(103.00919 + 57.021464, 5, 8, 2, 2, 1), //CamC is a modification of cysteine
    I(113.08406, 6, 11, 1, 1, 0),
    //L(113.08406, 6, 13, 1, 1, 0),
    N(114.04293, 4, 6, 2, 2, 0),
    D(115.02694, 4, 5, 1, 3, 0),
    R(156.10111, 6, 12, 4, 1, 0),
    Q(128.05858, 5, 8, 2, 2, 0),
    K(128.09496, 6, 12, 2, 1, 0),
    E(129.04259, 5, 7, 1, 3, 0),
    M(131.04049, 5, 9, 1, 1, 1),
    H(137.05891, 6, 7, 3, 1, 0),
    F(147.06841, 9, 9, 1, 1, 0),
    //U(150.953636),  //too rare
    Y(163.06333, 9, 9, 1, 2, 0),
    W(186.07931, 11, 10, 2, 1, 0);
    //O(237.147727); //too rare

    private double mass;

    private final int[] atomCount;

    public int[] getAtomCount() {
        return atomCount;
    }

    private Acid(double mass, int carbon, int hydrogen, int nitrogen, int oxygen, int sulfur) {
        this.mass = mass;
        atomCount = new int[]{carbon, hydrogen, nitrogen, oxygen, sulfur};
    }

    public double getMass() {
        return mass;
    }

    public boolean match(double[] limits) {
       return limits[0] < mass && limits[1] > mass;
    }

    public static HashMap<Character, Acid> acids = new HashMap<Character, Acid>();

    static {
        for (Acid acid : Acid.values()) {
            acids.put(acid.name().charAt(0), acid);
        }
    }

    public static Acid getAcid(char ch) {
        return acids.get(ch);
    }

    public static Acid getAcid(double mass) {
        double d = 1000000;
        Acid ans = null;
        for (Acid acid : acids.values()) {
            double newD = Math.abs(acid.getMass() - mass);
            if (newD < d) {
                d = newD;
                ans = acid;
            }
        }
        return ans;
    }

}
