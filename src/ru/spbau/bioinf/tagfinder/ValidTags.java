package ru.spbau.bioinf.tagfinder;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ValidTags {

    private Configuration conf;

    private int gap;

    public ValidTags(Configuration conf) {
        this.conf = conf;
    }

    private static NumberFormat df = NumberFormat.getInstance();
    static {
        df.setMaximumFractionDigits(2);
    }

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();
        Map<Integer,Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        ValidTags validTags = new ValidTags(conf);
        Set<Integer> usedProteins = new HashSet<Integer>();
        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            if (msAlignResults.containsKey(scanId)) {
                Integer proteinId = msAlignResults.get(scanId);
                if (usedProteins.contains(proteinId)) {
                    continue;
                }
                usedProteins.add(proteinId);
                if (scanId == 1946) {
                    //continue;
                }


                String sequence = proteins.get(proteinId).getSimplifiedAcids();
                String reverseSequence = getReverse(sequence);
                KDStatistics kdStatistics = new KDStatistics(conf);




                List<Peak> peaks =
                            scan.createStandardSpectrum()
                            //scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan))
                    ;


                //kdStatistics.generateEdges(peaks);
                //printUsualTagInfo(peaks, conf, scanId, proteinId, sequence, reverseSequence);

                validTags.gap = 3;
                kdStatistics.generateGapEdges(peaks, validTags.gap);
                validTags.printGappedTagInfo(peaks, scanId, proteinId, sequence, reverseSequence);
            }
        }
    }

    private static void printUsualTagInfo(List<Peak> peaks, Configuration conf, int scanId, Integer proteinId, String sequence, String reverseSequence) {
        Set<String> tags = new HashSet<String>();
        ValidTags validTags = new ValidTags(conf);
        validTags.generateTags(tags, peaks);

        Set<String> filteredTags = new HashSet<String>();
        for (String tag : tags) {
            String reverseTag = getReverse(tag);
            if (!filteredTags.contains(reverseTag) || reverseTag.equals(tag)) {
                filteredTags.add(tag);
            }
        }
        int[][] stat = new int[1000][2];
        for (String tag : filteredTags) {
            int len = tag.length();
            if (sequence.contains(tag) || reverseSequence.contains(tag)) {
                stat[len][0]++;
            } else {
                stat[len][1]++;
            }
        }
        printStat(stat, proteinId, scanId);
    }

    private void printGappedTagInfo(List<Peak> peaks, int scanId, Integer proteinId, String sequence, String reverseSequence) {
        int[][] stat = new int[1000][2];
        for (Peak peak : peaks) {
            Set<Integer> starts = new HashSet<Integer>();
            Set<Integer> reverseStarts = new HashSet<Integer>();
            for (int i = 0; i < sequence.length(); i++) {
                starts.add(i);
                reverseStarts.add(i);
            }
            processGappedTags(stat, peak, 0, sequence, reverseSequence, starts, reverseStarts);
        }
        printStat(stat, proteinId, scanId);
    }

    private void processGappedTags(int[][] stat, Peak peak, int prefix, String sequence, String reverseSequence, Set<Integer> starts, Set<Integer>reverseStarts) {
        if (prefix > 0) {
            stat[prefix][0]++;
        }
        List<Peak> nextPeaks = peak.getNext();
        for (Peak next : nextPeaks) {
            double[] limits = conf.getEdgeLimits(peak, next);
            Set<Integer> nextStarts = getNextStarts(sequence, starts, limits);
            Set<Integer> nextReverseStarts = getNextStarts(reverseSequence, reverseStarts, limits);
            if (nextStarts.size() + nextReverseStarts.size() == 0) {
                processWrongGappedTags(stat, next, prefix + 1);
            } else {
                processGappedTags(stat, next, prefix + 1, sequence, reverseSequence, nextStarts, nextReverseStarts);
            }
        }
    }

    private void processWrongGappedTags(int[][] stat, Peak peak, int prefix) {
        List<Peak> nextPeaks = peak.getNext();
        stat[prefix][1]++;
        for (Peak next : nextPeaks) {
            processWrongGappedTags(stat, next, prefix + 1);
        }
    }

    private Set<Integer> getNextStarts(String sequence, Set<Integer> starts, double[] limits) {
        Set<Integer> nextStarts = new HashSet<Integer>();
        for (int pos : starts) {
            if (pos >= sequence.length()) {
                continue;
            }
            double m = Acid.getAcid(sequence.charAt(pos)).getMass();
            if (limits[0] < m && m < limits[1]) {
                nextStarts.add(pos + 1);
            }
            if (pos  + 1 < sequence.length()) {
                m += Acid.getAcid(sequence.charAt(pos + 1)).getMass();
                if (limits[0] < m && m < limits[1]) {
                    nextStarts.add(pos + 2);
                }

            }
            if (gap > 2) {
                if (pos  + 2 < sequence.length()) {
                    m += Acid.getAcid(sequence.charAt(pos + 2)).getMass();
                    if (limits[0] < m && m < limits[1]) {
                        nextStarts.add(pos + 3);
                    }

                }
            }

        }
        return nextStarts;
    }

    private static void printStat(int[][] stat, Integer proteinId, int scanId) {
        System.out.print(scanId + " " +  proteinId);
        for (int i = 1; i < stat.length; i++) {
            int good = stat[i][0];
            int total = good + stat[i][1];
            if (total > 0) {
                System.out.print(" " + df.format((100d * good)/total));
            } else {
                break;
            }
        }
        System.out.println();
    }

    private static String getReverse(String tag) {
        return new StringBuilder(tag).reverse().toString();
    }


    public void generateTags(Set<String> tags, List<Peak> peaks) {
         for (Peak peak : peaks) {
             generateTags(tags, "", peak);
         }
     }

    public void generateTags(Set<String> tags, String prefix, Peak peak) {
        tags.add(prefix);
        for (Peak next : peak.getNext()) {
            for (Acid acid : Acid.values()) {
                if (acid.match(conf.getEdgeLimits(peak, next))){
                    generateTags(tags, prefix + acid.name(), next);
                }
            }
        }
    }


}