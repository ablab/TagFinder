package ru.spbau.bioinf.tagfinder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class CompareFilteredScans {

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        Configuration conf2 = new Configuration(args, UnmatchedScansGenerator.SHARED_MODE);

        List<Protein> proteins = conf.getProteins();
        Map<Integer,Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        Map<Integer, Scan> scans2 = conf2.getScans();
        KDStatistics kdStat = new KDStatistics(conf);
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Set<Integer> usedProteins = new HashSet<Integer>();
        Collections.sort(keys);
        for (int key : keys) {
            Scan scan2 = scans2.get(key);
            if (scan2 != null) {
                Scan scan = scans.get(key);
                int scanId = scan.getId();
                if (msAlignResults.containsKey(scanId)) {
                    Integer proteinId = msAlignResults.get(scanId);
                    if (usedProteins.contains(proteinId)) {
                        continue;
                    }
                    usedProteins.add(proteinId);
                    String sequence = proteins.get(proteinId).getSimplifiedAcids();
                    KD kd = kdStat.findKdBetweenBY(scan, sequence);
                    List<Peak> leftPeaks = scan2.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan));
                    GraphUtil.generateEdges(conf, leftPeaks);
                    Peak[] bestTag = GraphUtil.findBestTag(leftPeaks);
                    if (bestTag.length > 3) {
                        System.out.println(scanId + " " + proteinId + " " +  " " +  kd.toString() + " " + bestTag.length);
                    }
                }
            }
        }
    }
}
