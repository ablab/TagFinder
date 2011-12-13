package ru.spbau.bioinf.tagfinder;


import edu.ucsd.msalign.align.prsm.PrSM;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
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
    private static Map<Integer, Integer> ans = new HashMap<Integer, Integer>();

    public static void main(String[] args) throws Exception {
        File cacheFile = new File("cache.txt");
        BufferedReader in = ReaderUtil.createInputReader(cacheFile);
        cacheOut = new PrintWriter(new OutputStreamWriter(
                new FileOutputStream(cacheFile, true), "UTF-8"));
        do {
            String  s = in.readLine();
            if (s == null) {
                break;
            }
            String[] data = s.split(" ");
            evalues.put(data[0] + "_" + data[1], Double.parseDouble(data[2]));
        } while (true);

        int count = 0;
        Configuration conf = new Configuration(args);
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        System.out.println("MS-Align results: " + msAlignResults.keySet().size());
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
            double eValue = getEValueWrapper(scanId, proteinId);
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

        List<Integer> unmatchedScans2 = new ArrayList<Integer>();
        for (int scanId : unmatchedScans) {
            if (checkScanAgainstProteins(scanId, discoveredProteins)) {
                count++;
            } else {
                unmatchedScans2.add(scanId);
            }
        }

        long finish2 = System.currentTimeMillis();
        System.out.println(count + " matches  for plus " + (finish2 - finish));


        List<Integer> unmatchedScans3 = new ArrayList<Integer>();
        for (int scanId : unmatchedScans2) {
            Scan scan = scans.get(scanId);
            List<Peak> spectrum = scan.createSpectrumWithYPeaks(0);
            GraphUtil.generateEdges(conf, spectrum);
            Set<String> tags = GraphUtil.generateTags(conf, spectrum);
            List<Integer> proteinsForCheck = new ArrayList<Integer>();
            for (int len = 10; len >=5; len--) {
                for (String tag : tags) {
                    if (tag.length() == len) {
                        for (Protein protein : proteins) {
                            if (protein.getSimplifiedAcids().contains(tag)) {
                                int proteinId = protein.getProteinId();
                                if (!proteinsForCheck.contains(proteinId)) {
                                       proteinsForCheck.add(proteinId);
                                       if (proteinsForCheck.size() > 100) {
                                           break;
                                       }
                                   }
                            }
                        }
                    }
                }
            }

            if (checkScanAgainstProteins(scanId, proteinsForCheck)) {
                count++;
            } else {
                unmatchedScans3.add(scanId);
            }
        }

        long finish3 = System.currentTimeMillis();
        System.out.println(count + " matches  for plus " + (finish3 - finish2));

        //cacheOut.println("end matches");
        cacheOut.close();
        System.out.println("PrSM found " + ans.keySet().size());
        Map<Integer, Double> evaluesOld = conf.getEvalues();
        System.out.println("MS-Align only");
        for (Integer scanId : msAlignResults.keySet()) {
            if (!ans.containsKey(scanId)) {
                int proteinId = msAlignResults.get(scanId);
                System.out.println(scanId + " " + proteinId + " " + evaluesOld.get(scanId) + " " + getEValue(scanId, proteinId));
            }
        }
        System.out.println("New matches ");
        for (Integer scanId : ans.keySet()) {
            if (!msAlignResults.containsKey(scanId)) {
                int proteinId = ans.get(scanId);
                System.out.println(scanId + " " + proteinId + " " + getEValue(scanId, proteinId));
            }
        }

    }

    private static  boolean checkScanAgainstProteins(int scanId, Collection<Integer> proteins) throws Exception {
        for (int proteinId : proteins) {
            double eValue = getEValueWrapper(scanId, proteinId);
            if (eValue < Configuration.EVALUE_LIMIT) {
                System.out.println(scanId + " " + proteinId + " " + eValue);
                return true;
            }
        }
        return false;
    }

    private static double getEValueWrapper(int scanId, int proteinId) throws Exception {
        double ret = getEValue(scanId, proteinId);
        if (ret < Configuration.EVALUE_LIMIT) {
            ans.put(scanId, proteinId);
        }
        return ret;
    }
    
    private static double getEValue(int scanId, int proteinId) throws Exception {
        String key = scanId + "_" + proteinId;
        double ans = Integer.MAX_VALUE;
        if (evalues.containsKey(key)) {
            return evalues.get(key);
        } else {
            PrSM prsm = EValueAdapter.getBestEValue(scans.get(scanId), proteinId);
            if (prsm != null) {
                ans = prsm.getEValue();
            }
            evalues.put(key, ans);
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
