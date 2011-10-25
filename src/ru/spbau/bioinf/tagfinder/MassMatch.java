package ru.spbau.bioinf.tagfinder;

public class MassMatch implements Comparable<MassMatch> {

    private double mass;
    private int begin;
    private int end;
    private double mod;
    private double error;
    private double relativeError;

    public MassMatch(double mass, int begin, int end, double mod, double[] proteinSpectrum) {
        this.mass = mass;
        this.begin = begin;
        this.end = end;
        this.mod = mod;
        double diff = proteinSpectrum[end] - proteinSpectrum[begin];
        double adjustMass = mass + mod;
        error = adjustMass - diff;
        relativeError = Math.abs(error/diff);
    }

    public int compareTo(MassMatch other) {
        if (begin < other.begin) {
            return -1;
        } else if (begin > other.begin) {
            return 1;
        } else {
            return other.getMass() - mass < 0 ? 1 : - 1;
        }
    }

    public double getMass() {
        return mass;
    }

    public int getBegin() {
        return begin;
    }

    public int getEnd() {
        return end;
    }

    public double getMod() {
        return mod;
    }

    public double getError() {
        return error;
    }

    public double getRelativeError() {
        return relativeError;
    }
}
