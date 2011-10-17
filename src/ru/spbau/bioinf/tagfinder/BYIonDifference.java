package ru.spbau.bioinf.tagfinder;


import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class BYIonDifference {

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
                double bestShiftB = ShiftEngine.getBestShift(peaks, proteinSpectrum);
                double scoreB = ShiftEngine.getScore(ShiftEngine.getSpectrum(peaks), proteinSpectrum, bestShiftB);
                List<Peak> yPeaks = scan.getYPeaks();
                double bestShiftY = ShiftEngine.getBestShift(yPeaks, proteinSpectrum);
                double scoreY = ShiftEngine.getScore(ShiftEngine.getSpectrum(yPeaks), proteinSpectrum, bestShiftY);
                System.out.println(scanId + " " + proteinId + " " + df.format(bestShiftB) + " " + scoreB + " " + df.format(bestShiftY) + " " + scoreY + " " + (df.format(bestShiftY -bestShiftB)));
            }
        }
    }
}
