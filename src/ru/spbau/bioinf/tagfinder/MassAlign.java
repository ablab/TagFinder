package ru.spbau.bioinf.tagfinder;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import ru.spbau.bioinf.tagfinder.stat.PlaceStatistics;

public class MassAlign {

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
        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            if (msAlignResults.containsKey(scanId)) {

                Integer proteinId = msAlignResults.get(scanId);
                Protein protein = proteins.get(proteinId);
                double[] proteinSpectrum = ShiftEngine.getSpectrum(protein.getSimplifiedAcids());
                List<Peak> peaks = scan.getPeaks();
                List<MassMatch> massAlign = getMassAlign(proteinSpectrum, peaks);
                double count = 0;
                for (MassMatch massMatch : massAlign) {
                    if (massMatch.getRelativeError() < PlaceStatistics.TEN_PPM) {
                        count++;
                    }
                }

                System.out.println(scanId + " " + proteinId + " " + (int)count + " " + peaks.size() + " " + (count/peaks.size()));

            }
        }
    }

    public static List<MassMatch> getMassAlign(double[] proteinSpectrum, List<Peak> peaks) {
        List<MassMatch> ans = new ArrayList<MassMatch>();
        double[] mods = new double[]{0, -Consts.WATER, - Consts.AMMONIA, -1, +1, Consts.AMMONIA, Consts.WATER};
        for (Peak peak : peaks) {
            double mass = peak.getMass();
            MassMatch match = new MassMatch(mass, 0, proteinSpectrum.length - 1, 0, proteinSpectrum);
            for (int i = 0; i < proteinSpectrum.length; i++) {
                for (int j = i+1; j < proteinSpectrum.length; j++) {
                    for (int k = 0; k < mods.length; k++) {
                        double mod = mods[k];
                        MassMatch nextMatch = new MassMatch(mass, i, j, mod, proteinSpectrum);
                        if (nextMatch.getRelativeError() < match.getRelativeError()) {
                            match = nextMatch;
                        }
                    }
                }
            }
            ans.add(match);
        }
        Collections.sort(ans);
        return ans;
    }
}
