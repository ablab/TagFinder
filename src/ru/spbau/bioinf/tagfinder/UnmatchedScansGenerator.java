package ru.spbau.bioinf.tagfinder;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class UnmatchedScansGenerator {
    private Configuration conf;

    public UnmatchedScansGenerator(Configuration conf) {
        this.conf = conf;
    }

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();
        Map<Integer,Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);


        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            if (msAlignResults.containsKey(scanId)) {
                int proteinId = msAlignResults.get(scanId);
                String sequence = proteins.get(proteinId).getSimplifiedAcids();
                String reverseSequence = ValidTags.getReverse(sequence);

                List<Peak> peaks = scan.createStandardSpectrum();
                List<Peak> removed = new ArrayList<Peak>();
                for (Peak p1 : peaks) {
                    for (Peak p2 : p1.getNext()) {
                        for (Acid a1 : Acid.values()) {
                            if (a1.match(conf.getEdgeLimits(p1, p2))) {
                                for (Peak p3 : p2.getNext()) {
                                    for (Acid a2 : Acid.values()) {
                                        if (a2.match(conf.getEdgeLimits(p2, p3))) {
                                            for (Peak p4 : p3.getNext()) {
                                                for (Acid a3 : Acid.values()) {
                                                    if (a3.match(conf.getEdgeLimits(p3, p4))) {
                                                        String tag = a1.name() + a2.name() + a3.name();
                                                        if (sequence.contains(tag) || reverseSequence.contains(tag)) {
                                                            removed.add(p1);
                                                            removed.add(p2);
                                                            removed.add(p3);
                                                            removed.add(p4);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                for (Peak peak : peaks) {
                    if (peak.getIntensity() == 0) {
                        removed.add(peak);
                    }
                }
                peaks.removeAll(removed);
                Scan filtered = new Scan(scan, peaks);
                filtered.save(new File(conf.getInputDir(), "env2"));
            }
        }
    }
}
