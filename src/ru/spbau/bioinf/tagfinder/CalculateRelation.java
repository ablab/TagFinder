package ru.spbau.bioinf.tagfinder;

import java.io.BufferedReader;
import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class CalculateRelation {

    public static final int MAX_TAG = 100;

    private static NumberFormat df = NumberFormat.getInstance();
    static {
        df.setMaximumFractionDigits(2);
    }

    public static void main(String[] args) throws Exception {
        for (int gap = 1; gap <=3; gap++) {
            compare("share_bar_basic_" + gap + "_correct.txt", "share_bar_basic_" + gap +"_proper.txt");
            correctD("share_bar_basic_" + gap +"_correct.txt");
            correctD("share_bar_virtual_full_" + gap +"_correct.txt");
        }
    }

    public static void compare(String fileFirst, String fileSecond) throws Exception {
        System.out.println(fileFirst + " " + fileSecond);
        BufferedReader inOne = ReaderUtil.createInputReader(new File("res", fileFirst));
        BufferedReader inTwo = ReaderUtil.createInputReader(new File("res", fileSecond));
        List<long[]> pairs = new ArrayList<long[]>();
        long[][][] stat = new long[10000][MAX_TAG][4];
        int n = 0;
        List<Double>[] goodToAll = new List[100];
        List<Double>[] goodToGood = new List[100];
        List<Double>[] allToAll = new List[100];
        for (int i = 0; i < MAX_TAG; i++) {
            goodToAll[i] = new ArrayList<Double>();
            goodToGood[i] = new ArrayList<Double>();
            allToAll[i] = new ArrayList<Double>();
        }

        do {
            String s1 = inOne.readLine();
            String s2 = inTwo.readLine();
            if (s1.indexOf(",") > 0) {
                break;
            }
            long[] d1 = getData(s1);
            long[] d2 = getData(s2);
            pairs.add(new long[] {d1[0], d1[1]});
            int pos = 2;
            int d = 1;

            while (pos < d1.length) {
                stat[n][d][0] = d1[pos];
                stat[n][d][1] = d1[pos + 1];
                d++;
                pos += 2;
            }

            d = 1;
            pos = 2;
            while (pos < d2.length) {
                stat[n][d][2] = d2[pos];
                stat[n][d][3] = d2[pos + 1];
                d++;
                pos += 2;
            }

            d = 1;
            do {
                long[] q = stat[n][d];
                double total =  q[2] + q[3];
                if (total == 0) {
                    break;
                }
                goodToAll[d].add((q[0])/total);
                allToAll[d].add((q[0] + q[1])/total);
                if (q[2] > 0) {
                    double v = q[2];
                    goodToGood[d].add(q[0]/v);
                }

                d++;
            } while (true);
            n++;
        } while(true);
        System.out.println("Good to All: ");
        printPercentage(goodToAll);
        System.out.println("Good to Good: ");
        printPercentage(goodToGood);
        System.out.println("All to All: ");
        printPercentage(allToAll);
    }

    public static void correctD(String file) throws Exception {
        System.out.println(file);
        BufferedReader in = ReaderUtil.createInputReader(new File("res", file));

        List<long[]> pairs = new ArrayList<long[]>();
        int n = 0;
        int[] dShare = new int[MAX_TAG];
        do {
            String s1 = in.readLine();

            if (s1.indexOf(",") > 0) {
                break;
            }
            long[] d1 = getData(s1);
            pairs.add(new long[] {d1[0], d1[1]});
            int pos = 2;
            int d = 1;

            int maxD = 0;
            while (pos < d1.length) {
                if (d1[pos] > 0) {
                    maxD = d;
                }
                d++;
                pos += 2;
            }
            dShare[maxD]++;
        } while(true);
        double total = pairs.size();
        System.out.println("Distribution of longest tags: ");
        for (int i = 0; i < dShare.length; i++) {
            int v = dShare[i];
            if (v > 0) {
                System.out.println(i + " " + v + " " + df.format(100 * v / total));
            }

        }

    }

    private static void printPercentage(List<Double>[] stat) {
        for (int d = 1; d < MAX_TAG; d++) {
            List<Double> values = stat[d];
            if (values.size() == 0) {
                break;
            }
            double total = 0;
            for (double value : values) {
                total += value;
            }
            System.out.print(df.format(100 * total / values.size()) + " ");
        }
        System.out.println();
    }

    public static long[] getData(String s) {
        String[] d = s.split(" ");
        long[] ans = new long[d.length];
        for (int i = 0; i < ans.length; i++) {
            ans[i] = Long.parseLong(d[i]);
        }
        return ans;
    }
}
