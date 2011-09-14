package ru.spbau.bioinf.tagfinder;

public class KD implements Comparable<KD> {
    private int k;
    private int d;

    public KD(int k, int d) {
        this.k = k;
        this.d = d;
    }

    public int getK() {
        return k;
    }

    public int getD() {
        return d;
    }

    @Override
    public int hashCode() {
        return k*1024 + d;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof KD)) {
            return false;
        }
        KD kd = (KD) obj;
        return k == kd.getK() && d == kd.getD();
    }

    @Override
    public String toString() {
        return "(" + k + ", " + d + ")";
    }

    @Override
    public int compareTo(KD other) {
        if (other.getK() > k) {
            return 1;
        } else if (other.getK() < k) {
            return -1;
        } else if (other.getD() > d) {
            return 1;
        } else if (other.getD() < d) {
            return -1;
        }
        return 0;
    }
}
