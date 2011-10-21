package ru.spbau.bioinf.tagfinder;


import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class ShiftsDifference {

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
                List<Double> shifts = ShiftEngine.getBestShifts(peaks, scan.getPrecursorMass(), proteinSpectrum);
                if (shifts.size() > 1) {
                    System.out.print(scanId + " " + proteinId);
                    double[] spectrum = ShiftEngine.getSpectrum(scan.createSpectrumWithYPeaks(scan.getPrecursorMass()));
                    for (Double shift : shifts) {
                        double score = ShiftEngine.getScore(spectrum, proteinSpectrum, shift);
                        System.out.print(" " + df.format(shift) + " " + (int) score);
                    }
                    System.out.println();
                    for (int i = 0; i < shifts.size(); i++) {
                        for (int j = i + 1; j < shifts.size(); j++) {
                            double delta = shifts.get(j) - shifts.get(i);
                            for (Peak peak : peaks) {
                                if (Math.abs(Math.abs(delta) - peak.getMass()) <0.1) {
                                    System.out.println(shifts.get(i) + " " + shifts.get(j) + " " +  peak.getMass() + " " + peak.getIntensity() + " " + peak.getCharge());
                                }
                            }

                        }
                    }
                }
            }
        }
    }
}
