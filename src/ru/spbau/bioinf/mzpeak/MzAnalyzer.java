package ru.spbau.bioinf.mzpeak;

import java.io.File;
import java.util.List;
import java.util.Map;
import ru.spbau.bioinf.pfind.Finder;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.Consts;
import ru.spbau.bioinf.tagfinder.Protein;

public class MzAnalyzer {

    public static void main(String[] args) throws Exception {
        Map<Integer,MzScan> mzScans = MzReader.getMzScans();
        System.out.println(mzScans.size());
        Configuration conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();

        /*
        MzScan scan = null;
        for (MzScan mzScan : scans) {
            if (mzScan.getScanId() == 3456) {
                scan = mzScan;
                break;
            }
        }
        System.out.println(scan.getScanId() + " " + scan.getPoints().size());
        */
        
        //plusMinusCheck(scans);
        //for (int i = 1; i < 10; i++) {
        int i = 1;

        int count = 0;
        for (MzScan scan : mzScans.values()) {
            List<Double> points = scan.getSupportedPoints(0.0000001d * i, 0);

            int bestScore = 0;
            int proteinId = 0;
            for (Protein protein : proteins) {
                int score = getScore(points, protein);
                if (score > bestScore) {
                    bestScore = score;
                    proteinId = protein.getProteinId();
                }
            }
            if (bestScore > 0) {
                try {
                    double eValue = Finder.getEValue(scan.getScanId(), proteinId);
                    if (eValue < 0.0024) {
                        count++;
                    }
                    System.out.println(scan.getScanId() + " " +  proteinId + " " + eValue + " " + bestScore + " " + count);
                } catch (Throwable e) {
                }
            }
        }

        System.out.println("count = " + count);

        //System.out.println("scan.getSupportedPoints " + i + " " + points.size());
        //}
        for (MzScan mzScan : mzScans.values()) {
            //getCommonPeaks(scan, mzScan);
        }
    }

    private static int getScore(List<Double> points, Protein protein) {
        double[] suffixes = protein.getYEnds();
        int score = 0;
        for (double suffix : suffixes) {                
            for (double point : points) {
                if (Math.abs(point - suffix - Consts.WATER) < 0.1) {
                    score++;
                }
            }
        }
        return score;
    }

    private static void plusMinusCheck(List<MzScan> scans) {
        int t1 = 0;
        int t2 = 0;
        int t3 = 0;
        for (MzScan mzScan : scans) {
            int s1 = mzScan.getSupportedPoints(0.0000001, -1).size();
            int s2 = mzScan.getSupportedPoints(0.0000001, +1).size();

            if (s1 + s2 > 10) {
                if (s1 > s2) {
                    t1++;
                } else if (s2 > s1){
                    t2++;
                } else {
                    t3++;
                }

                System.out.println(mzScan.getScanId() + " " + s1 + " " + s2);
            }
        }

        System.out.println("t1 = " + t1);
        System.out.println("t2 = " + t2);
        System.out.println("t3 = " + t3);
    }

    private static void getCommonPeaks(MzScan scan, MzScan mzScan) {
        List<MzPoint> p2 = mzScan.getPoints();
        List<MzPoint> p1 = scan.getPoints();
        int cur1 = 0;
        int cur2 = 0;
        int common = 0;
        while (cur1 < p1.size() && cur2 < p2.size()) {
            double m1 = p1.get(cur1).getMass();
            double m2 = p2.get(cur2).getMass();
            if (Math.abs(m1-m2) < 0.000002 * (m1 + m2)) {
                cur1++;
                cur2++;
                common++;
            } else if (m1 < m2) {
                cur1++;
            } else {
                cur2++;
            }
        }
        if (common > 500) {
            System.out.println(mzScan.getScanId() + " " + common + " " + (p2.size() - common));
        }
    }
}
