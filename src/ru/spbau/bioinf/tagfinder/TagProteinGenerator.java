package ru.spbau.bioinf.tagfinder;

import edu.ucsd.msalign.align.prsm.PrSM;
import java.util.HashMap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class TagProteinGenerator {

    private static Configuration conf = new Configuration(new String[]{}, UnmatchedScansGenerator.SHARED_MODE);

    private static Map<Integer,Integer> msAlignResults;


    public static void main(String[] args) throws Exception {
        List<Protein> proteins = conf.getProteins();
        Map<Integer, Integer> badMSAlignResults = conf.getBadMSAlignResults();
        msAlignResults = new Configuration(args).getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        List<String> sequences = new ArrayList<String>();
        for (Protein protein : proteins) {
            sequences.add(protein.getSimplifiedAcids());
        }

        EValueAdapter.init(conf);


        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            int oldProteinId = -1;

            if (msAlignResults.containsKey(scanId)) {
                //continue;
                oldProteinId = msAlignResults.get(scanId);
            }

            if (badMSAlignResults.containsKey(scanId)) {
                oldProteinId = badMSAlignResults.get(scanId);
            }
            if (oldProteinId == -1) {
                //continue;
            }

            //System.out.println("matchedProteinId = " + matchedProteinId);

            List<Peak> peaks = //scan.createStandardSpectrum();
                    scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan));
            GraphUtil.generateEdges(conf, peaks);
            generateFiveAcidsTags(peaks, sequences, oldProteinId, scan);
        }
    }

    private static Set<String> done = new HashSet<String>();

    public static void generateFiveAcidsTags(List<Peak> peaks, List<String> sequences, int oldProteinId, Scan scan) throws Exception {
        List<Set<String>> allTags = new ArrayList<Set<String>>();
        for (int depth = 10; depth >= 5; depth--) {
            Set<String> tags = getTags(peaks, depth);
            if (tags.size() > 0) {
                allTags.add(tags);
            }
        }

        HashMap<String, Integer> stat = new HashMap<String, Integer>();
        for (Set<String> tagSet : allTags) {
            for (String tag : tagSet) {
                int count = 0;
                for (int i = 0; i < sequences.size(); i++) {
                    String seq = sequences.get(i);
                    if (seq.contains(tag)) {
                        count++;
                    }
                }
                stat.put(tag, count);
            }
        }



        for (int proteinId = 0; proteinId < sequences.size(); proteinId++) {
            String seq = sequences.get(proteinId);
            Set<String> found = new HashSet<String>();
            boolean isMatch = false;
            for (Set<String> tagSet : allTags) {
                for (String tag : tagSet) {
                    boolean isSubTag = false;
                    for (String oldTag : found) {
                        if (oldTag.contains(tag)) {
                            isSubTag = true;
                            break;
                        }
                    }
                    if (isSubTag) {
                        continue;
                    }
                    if (seq.contains(tag)) {
                        found.add(tag);
                        if (found.size() > 2 || tag.length() >= 6) {
                            isMatch = true;
                        }
                    }
                }
            }
            if (isMatch) {
                Set<String> bestTags = allTags.get(0);
                int maxTag = bestTags.iterator().next().length();
                String scanText = "Scan " + scan.getId() + " with " + bestTags.size() + " " + maxTag + "-aa tags ";
                if (oldProteinId >= 0) {
                    scanText += " matched to protein " + oldProteinId;
                }

                int otherMatches = 0;
                for (Map.Entry<Integer, Integer> entry : msAlignResults.entrySet()) {
                    if (entry.getValue().equals(proteinId)) {
                        otherMatches++;
                    }
                }

                if (otherMatches == 0) {
                    //continue;
                }

                String key = scan.getId() + "_" + proteinId;


                if (!done.contains(key)) {
                    PrSM prsm = EValueAdapter.getBestEValue(scan, proteinId);
                    double eValue = prsm == null? 9E10 : prsm.getEValue();
                    done.add(key);
                    System.out.println(scan.getId() + " " + oldProteinId + " " + proteinId + " " + eValue);

                }


            }
        }

    }

    public static Set<String> getTags(List<Peak> peaks, int depth) {
        Set<String> tags = new HashSet<String>();
        for (Peak peak : peaks) {
            generateTags(tags, peak, "", depth);
        }
        return tags;
    }

    public static void generateTags(Set<String> used, Peak peak, String prefix, int depth) {
        //System.out.println(prefix + " " + peak.getValue() + " " + peak.getMass());
        if (prefix.length() == depth) {
            used.add(prefix);
        } else {
            for (Peak next : peak.getNext()) {
                for (Acid acid : Acid.values()) {
                    if (acid.match(conf.getEdgeLimits(peak, next))) {
                        generateTags(used, next, prefix + acid.name(), depth);
                    }
                }
            }

        }
    }

}
