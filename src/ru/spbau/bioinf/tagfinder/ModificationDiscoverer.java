package ru.spbau.bioinf.tagfinder;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import ru.spbau.bioinf.tagfinder.stat.PlaceStatistics;

public class ModificationDiscoverer {

    private static NumberFormat df = NumberFormat.getInstance();

    static {
        df.setGroupingUsed(false);
        df.setMaximumFractionDigits(2);
    }

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);

        //keys.clear();keys.add(556);

        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            if (msAlignResults.containsKey(scanId)) {
                Integer proteinId = msAlignResults.get(scanId);
                Protein protein = proteins.get(proteinId);
                double[] proteinSpectrum = ShiftEngine.getSpectrum(protein.getSimplifiedAcids());
                List<Peak> peaks = scan.getPeaks();

                getBestDiff(scanId, proteinId, proteinSpectrum,
                        getBadMatches(MassAlign.getMassAlign(proteinSpectrum, peaks)));
            }
        }
    }

    public static List<MassMatch> getBadMatches(List<MassMatch> matches) {
        List<MassMatch> ans = new ArrayList<MassMatch>();
        for (MassMatch match : matches) {
            if (match.getRelativeError() > PlaceStatistics.TEN_PPM) {
                ans.add(match);
            }
        }
        return ans;
    }

    public static double getBestDiff(int scanId, Integer proteinId, double[] proteinSpectrum, List<MassMatch> badMatches) {
        List<Double> diffs = new ArrayList<Double>();
        for (int i = 0; i < proteinSpectrum.length; i++) {
            for (int j = i+1; j < proteinSpectrum.length; j++) {
                diffs.add(proteinSpectrum[j] - proteinSpectrum[i]);
            }
        }

        double[] d = new double[diffs.size()];
        for (int i = 0; i < d.length; i++) {
            d[i] = diffs.get(i);
        }
        d = ShiftEngine.merge(d);

        //Collections.sort(diffs);
        List<Double> allDiffs = new ArrayList<Double>();

        for (MassMatch massMatch : badMatches) {
            for (Double diff : d) {
                allDiffs.add(diff - massMatch.getMass());
            }
        }

        Collections.sort(allDiffs);
        if (allDiffs.size() == 0) {
            return 0;
        }
        double prev = allDiffs.get(0) - 100;
        int count = 0;
        double min = prev;
        int bestScore = 0;
        double bestValue = 0;

        for (double s : allDiffs) {
            if (s - prev < 0.01) {
                count++;
                prev = s;
            } else {
                if (count > 10) {
                    double nextBestValue = (prev + min) / 2;
                    if (Math.abs(nextBestValue) < 40) {
                        if (bestScore < count) {
                            bestScore = count;
                            bestValue = nextBestValue;
                            System.out.println(scanId + " " + proteinId + " bestScore = " + bestScore + " " + bestValue);
                        }
                    }
                }
                count = 1;
                prev = s;
                min = s;
            }
        }
        return bestValue;
    }

}
