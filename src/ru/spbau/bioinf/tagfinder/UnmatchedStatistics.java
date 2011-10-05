package ru.spbau.bioinf.tagfinder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class UnmatchedStatistics {

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        Map<Integer,Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        int[] stat = new int[100];
        int max = 0;
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            if (!msAlignResults.containsKey(scanId)) {
                List<Peak> peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan));
                GraphUtil.generateEdges(conf, peaks);
                Peak[] bestTag = GraphUtil.findBestTag(peaks);
                int v = bestTag.length;
                System.out.println(scanId + " " + v);
                stat[v]++;
                if (max < v) {
                    max = v;
                }
            }
        }
        for (int i = 0; i <= max; i++) {
            System.out.println(i + " " + stat[i]);
        }
    }
}
