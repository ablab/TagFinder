package ru.spbau.bioinf.tagfinder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ShiftCheck {

    public static final String GOOD = "good";
    public static final String OTHER = "other";
    private static List<Protein> proteins;

    public static void main(String[] args) throws Exception {

        Configuration conf = new Configuration(args);

        proteins = conf.getProteins();
        Map<Integer, Scan> scans = conf.getScans();
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(msAlignResults.keySet());
        Collections.sort(keys);
        int good = 0;
        int bad = 0;
        int yellow = 0;

        Map<Integer, Double> evalues = conf.getEvalues();
        Set<Integer> matchedProteins = new HashSet<Integer>();
        Set<Integer> unmatchedScans = new HashSet<Integer>();

        for (int scanId : keys) {
            Scan scan = scans.get(scanId);
            List<Peak> peaks = scan.getPeaks();
            Collections.sort(peaks);
            double[] p = new double[peaks.size()];
            for (int i = 0; i < p.length; i++) {
                p[i] = peaks.get(i).getMass(); //+ shift;
            }
            int proteinId = msAlignResults.get(scanId);
            Protein protein = proteins.get(proteinId);
            double[] yEnds = protein.getYEnds();
            double[] bEnds = protein.getBEnds();
            if (yEnds.length < 1 || p.length < 1) {
                continue;
            }
            //int score = scoring(p, yEnds, bEnds);
            //if (score < 27) {

                int newProteinId = getBestProtein(p);
                if (newProteinId != proteinId) {
                    unmatchedScans.add(scanId);
                    bad++;
                } else {
                    matchedProteins.add(proteinId);
                    yellow++;
                }
            //} else {
            //    good++;
            //    matchedProteins.add(proteinId);
            //}
        }
        int old = 0;
        int newBad = 0;
        for (int scanId : unmatchedScans) {
            int proteinId = msAlignResults.get(scanId);
            if (matchedProteins.contains(proteinId)) {
                old++;
            } else {
                System.out.println(scanId + " " + proteinId + " " + scans.get(scanId).getPeaks().size() + " " + evalues.get(scanId));
                newBad++;
            }
        }

        System.out.println("total = " + scans.keySet().size());
        System.out.println("matched = " + msAlignResults.keySet().size());
        //System.out.println("good = " + good);
        System.out.println("yellow = " + yellow);
        System.out.println("bad = " + bad);
        System.out.println("old = " + old);
        System.out.println("newBad = " + newBad);




        for (int scoreLimit = 35; scoreLimit >= 25; scoreLimit--) {
            //process(scoreLimit, conf);
            //System.out.println("============================");
        }

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
            int score = scoring(p, yEnds, bEnds);
                    /*
                    max(
                    getScoreY(p, yEnds),
                    getScoreB(p, bEnds),
                    getScoreBM(p, bEnds)

            );        */

            if (score > bestScore) {
                bestScore = score;
                ans = protein.getProteinId();
            }
        }
        return ans;
    }

    private static int scoring(double[] p, double[] yEnds, double[] bEnds) {
        return max(
                getScoreY(p, yEnds),
                getScoreB(p, bEnds),
                getScoreBM(p, bEnds),
                getScoreBN(p, bEnds)

        );
    }

    private static void process(int scoreLimit, Configuration conf) throws Exception {
        List<Protein> proteins = conf.getProteins();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();

        System.out.println("scoreLimit = " + scoreLimit);
        int good = 0;
        int bad = 0;
        int other = 0;
        EValueAdapter.init(conf);
        int totalFail = 0;
        int onlyGood = 0;
        int badGoodScore = 0;
        for (int scanId : keys) {
            Scan scan = scans.get(scanId);
            StringBuilder text = new StringBuilder();
            List<Peak> peaks = scan.getPeaks();
            Collections.sort(peaks);
            double[] p = new double[peaks.size()];
            for (int i = 0; i < p.length; i++) {
                p[i] = peaks.get(i).getMass(); //+ shift;
            }

            boolean noGood = true;
            boolean hasFail = false;
            boolean printResult = false;
            int goodScore = 0;
            int maxBadScore = 0;
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
                        getScoreBM(p, bEnds)

                );

                if (score >= scoreLimit) {
                    int proteinId = protein.getProteinId();
                    String res = "other";
                    int oldProteinId = 0;
                    if (msAlignResults.containsKey(scanId)) {
                        oldProteinId = msAlignResults.get(scanId);
                        res = oldProteinId == proteinId ? GOOD : "fail " + oldProteinId;
                    }

                    text.append(scanId + " " + proteinId + " " + score + " " + res + " ");

                    if (res.startsWith(OTHER)) {
                        other++;
                        printResult = true;
                        text.append(EValueAdapter.getBestEValue(scan, proteinId).getEValue());
                    } else if (res.startsWith(GOOD)) {
                        good++;
                        System.out.println(scanId + " " + proteinId + " " + score);
                        noGood = false;
                        goodScore = score;
                    } else {
                        bad++;
                        maxBadScore = Math.max(score, maxBadScore);
                        hasFail = true;
                        printResult = true;
                        text.append(EValueAdapter.getBestEValue(scan, proteinId).getEValue() + " " + EValueAdapter.getBestEValue(scan, oldProteinId).getEValue());
                    }
                    text.append("\n");
                }
            }
            if (noGood && hasFail) {
                System.out.println("TotalFail");
                totalFail++;
            }
            if (!noGood && !hasFail) {
                onlyGood++;
            }
            if (maxBadScore > 0 && maxBadScore >= goodScore) {
                text.append("BadGoodScore");
                badGoodScore++;
            }
            if (printResult) {
                System.out.println(text.toString());
            }

        }
        System.out.println("good = " + good);
        System.out.println("bad = " + bad);
        System.out.println("other = " + other);
        System.out.println("onlyGood = " + onlyGood);
        System.out.println("totalFail = " + totalFail);
        System.out.println("badGoodScore = " + badGoodScore);
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
