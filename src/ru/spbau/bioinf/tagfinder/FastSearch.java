package ru.spbau.bioinf.tagfinder;


import edu.ucsd.msalign.align.prsm.PrSM;
import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class FastSearch {

    private static List<Protein> proteins;
    private static Map<Integer,Scan> scans;

    private static Map<String, Double> evalues = new HashMap<String, Double>();
    private static PrintWriter cacheOut;

    public static void main(String[] args) throws Exception {
        BufferedReader in = ReaderUtil.createInputReader(new File("fs", "1.txt"));
        cacheOut = ReaderUtil.createOutputFile(new File("cacheOut.txt"));
        do {
            String  s = in.readLine();
            if (s.contains("matches")) {
                break;
            }
            String[] data = s.split(" ");
            evalues.put(data[0] + "_" + data[1], Double.parseDouble(data[2]));
        } while (true);

        int count = 0;
        Configuration conf = new Configuration(args);
        proteins = conf.getProteins();
        scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        EValueAdapter.init(conf);

        long start = System.currentTimeMillis();
        List<Integer> unmatchedScans = new ArrayList<Integer>();
        Set<Integer> discoveredProteins = new HashSet<Integer>();

        for (int scanId : keys) {
            Scan scan = scans.get(scanId);
            List<Peak> peaks = scan.getPeaks();
            Collections.sort(peaks);
            double[] p = new double[peaks.size()];
            for (int i = 0; i < p.length; i++) {
                p[i] = peaks.get(i).getMass(); //+ shift;
            }
            int proteinId = getBestProtein(p);
            double eValue = getEValue(scanId, proteinId);
            if (eValue < Configuration.EVALUE_LIMIT) {
                System.out.println(scanId + " " + proteinId + " " + eValue);
                count++;
                discoveredProteins.add(proteinId);
            } else {
                unmatchedScans.add(scanId);
            }
        }
        long finish = System.currentTimeMillis();
        System.out.println(count + " matches  found in " + (finish - start));

        for (Integer scanId : unmatchedScans) {
            for (int proteinId : discoveredProteins) {
                double eValue = getEValue(scanId, proteinId);
                if (eValue < Configuration.EVALUE_LIMIT) {
                    System.out.println(scanId + " " + proteinId + " " + eValue);
                    count++;
                }
            }
        }

        long finish2 = System.currentTimeMillis();
        System.out.println(count + " matches  for plus " + (finish2 - finish));

        cacheOut.println("end matches");
        cacheOut.close();

    }

    private static double getEValue(int scanId, int proteinId) throws Exception {
        String key = scanId + "_" + proteinId;
        double ans = Integer.MAX_VALUE;
        if (evalues.containsKey(key)) {
            ans = evalues.get(key);
        } else {
            PrSM prsm = EValueAdapter.getBestEValue(scans.get(scanId), proteinId);
            if (prsm != null) {
                ans = prsm.getEValue();
            }
        }
        cacheOut.println(scanId + " " + proteinId + " " + ans);
        cacheOut.flush();
        return ans;
    }

    public static int getBestProtein(double[] p) {
        int ans = 0;
        int bestScore = 0;
        for (Protein protein : proteins) {
            if (protein.getProteinId() == 2535) {
                continue;
            }
            double[] yEnds = protein.getYEnds();
            double[] bEnds = protein.getBEnds();
            if (yEnds.length < 1 || p.length < 1) {
                continue;
            }
            int score = max(
                    getScoreY(p, yEnds),
                    getScoreB(p, bEnds),
                    getScoreBM(p, bEnds),
                    getScoreBN(p, bEnds)
            );

            if (score > bestScore) {
                bestScore = score;
                ans = protein.getProteinId();
            }
        }
        return ans;
    }





    private static int max(int ... scores) {
        int ans = 0;
        for (int score : scores) {
            if (score > ans) {
                ans = score;
            }
        }
        return ans;
    }
    private static int getScoreY(double[] p, double[] yEnds) {
        return (int) Math.round(3 * ShiftEngine.getScore(p, yEnds, -Consts.WATER) + ShiftEngine.getScore(p, yEnds, 0));
    }

    private static int getScoreB(double[] p, double[] bEnds) {
        return (int) Math.round(ShiftEngine.getScore(p, bEnds, -Consts.WATER) + 3 * ShiftEngine.getScore(p, bEnds, 0));
    }

    private static int getScoreBM(double[] p, double[] bEnds) {
        return (int) Math.round(3 * ShiftEngine.getScore(p, bEnds, -Acid.M.getMass()) + ShiftEngine.getScore(p, bEnds, - Acid.M.getMass() - Consts.AMMONIA));
    }

    private static int getScoreBN(double[] p, double[] bEnds) {
        return (int) Math.round(3 * ShiftEngine.getScore(p, bEnds, -Consts.N));
    }

}
