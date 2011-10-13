package ru.spbau.bioinf.tagfinder;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ValidTags {

    private static final int MAX_GAPPED_TAG = 90;

    private Configuration conf;

    private int gap;
    private static Map<Integer,List<Peak>> msAlignPeaks;

    public ValidTags(Configuration conf) {
        this.conf = conf;
    }

    private static NumberFormat df = NumberFormat.getInstance();
    static {
        df.setMaximumFractionDigits(2);
    }

    private static boolean needReverseTag = false;

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
        msAlignPeaks = conf.getMSAlignPeaks();

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
                    continue;
                }


                Protein protein = proteins.get(proteinId);
                validTags.process(scan, protein, 3);
            }
        }
    }

    public void process(Scan scan, Protein protein, int gap) {
        this.gap = gap;
        String sequence = protein.getSimplifiedAcids();
        List<Peak> peaks =
                    //msAlignPeaks.get(scanId)
                    scan.createStandardSpectrum()
                    //scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan))
                    //scan.createStandardSpectrumWithOnes()
            ;

        needReverseTag = false;
        System.out.println("needReverseTag = " + needReverseTag);

        peaks = addOnes(peaks);

        //filterMonotags(peaks);
        int scanId = scan.getId();
        int proteinId = protein.getProteinId();

        if (gap == 1) {
            GraphUtil.generateEdges(conf, peaks);
            double[] positions = ShiftEngine.getPositions(peaks);
            printUsualTagInfo(peaks, conf, scanId, proteinId, sequence, positions);
        } else {
            GraphUtil.generateGapEdges(conf, peaks, gap);
            printGappedTagInfo(peaks, scanId, proteinId, sequence, getReverse(sequence));
        }
    }

    private List<Peak> addOnes(List<Peak> peaks) {
        List<Peak> ans = new ArrayList<Peak>();
        for (Peak peak : peaks) {
            ans.add(new Peak(peak.getValue() - 1, 0 ,0));
            ans.add(new Peak(peak.getValue(), 0 ,0));
            ans.add(new Peak(peak.getValue() + 1, 0 ,0));
        }
        Collections.sort(ans);
        return ans;
    }

    private static void filterMonotags(List<Peak> peaks) {
        for (Peak peak : peaks) {
            for (Iterator<Peak> iterator = peak.getNext().iterator(); iterator.hasNext(); ) {
                Peak next = iterator.next();
                if (next.getPeakType() != peak.getPeakType()) {
                    iterator.remove();
                }
            }
        }
    }

    private static void printUsualTagInfo(List<Peak> peaks, Configuration conf, int scanId, Integer proteinId, String sequence, double[] positions) {

        Set<String> tags = GraphUtil.generateTags(conf, peaks);

        long[][] stat = new long[1000][2];
        for (String tag : tags) {
            updateStat(stat, tag, sequence, peaks, positions);
        }
        printStat(stat, proteinId, scanId);
    }

    private static void updateStat(long[][] stat, String tag, String sequence, List<Peak> peaks, double[] positions) {
        String reverseTag = getReverse(tag);
        boolean updateForReverse = !reverseTag.equals(tag);
        for (double pos : positions) {
            if (GraphUtil.tagStartsAtPos(pos, tag, peaks)) {
                updateStat(stat, tag, sequence);
                if (updateForReverse) {
                    updateStat(stat, reverseTag, sequence);
                }
            }
        }
    }

    private static void updateStat(long[][] stat, String tag, String sequence) {
        int len = tag.length();
        if (sequence.contains(tag)) {
            stat[len][0]++;
        } else {
            stat[len][1]++;
        }
    }

    private void printGappedTagInfo(List<Peak> peaks, int scanId, Integer proteinId, String sequence, String reverseSequence) {
        long[][] stat = new long[1000][2];
        for (Peak peak : peaks) {
            Set<Integer> starts = new HashSet<Integer>();
            Set<Integer> reverseStarts = new HashSet<Integer>();
            for (int i = 0; i < sequence.length(); i++) {
                starts.add(i);
                if (needReverseTag) {
                    reverseStarts.add(i);
                }
            }
            processGappedTags(stat, peak, 0, sequence, reverseSequence, starts, reverseStarts);
        }
        printStat(stat, proteinId, scanId);
    }

    private void processGappedTags(long[][] stat, Peak peak, int prefix, String sequence, String reverseSequence, Set<Integer> starts, Set<Integer>reverseStarts) {
        if (prefix > 0) {
            stat[prefix][0]++;
        }
        if (prefix > MAX_GAPPED_TAG) {
            return;
        }

        List<Peak> nextPeaks = peak.getNext();
        for (Peak next : nextPeaks) {
            double[] limits = conf.getEdgeLimits(peak, next);
            Set<Integer> nextStarts = getNextStarts(sequence, starts, limits, gap);
            Set<Integer> nextReverseStarts = getNextStarts(reverseSequence, reverseStarts, limits, gap);
            if (nextStarts.size() + nextReverseStarts.size() == 0) {
                processWrongGappedTags(stat, next, prefix + 1);
            } else {
                processGappedTags(stat, next, prefix + 1, sequence, reverseSequence, nextStarts, nextReverseStarts);
            }
        }
    }

    private void processWrongGappedTags(long[][] stat, Peak peak, int prefix) {
        List<Peak> nextPeaks = peak.getNext();
        stat[prefix][1]++;
        long[] statOrig = new long[stat.length];
        for (int i = 0; i < statOrig.length; i++) {
            statOrig[i] = stat[i][1];
        }

        if (prefix > MAX_GAPPED_TAG) {
            return;
        }
        for (Peak next : nextPeaks) {
            processWrongGappedTags(stat, next, prefix + 1);
        }
    }

    public static Set<Integer> getNextStarts(String sequence, Set<Integer> starts, double[] limits, int gap) {
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

    private static void printStat(long[][] stat, Integer proteinId, int scanId) {
        System.out.print(scanId + " " +  proteinId);
        for (int i = 1; i < stat.length; i++) {
            long good = stat[i][0];
            long total = good + stat[i][1];
            if (total > 0) {
                System.out.print(" " + df.format((100d * good)/total));
            } else {
                break;
            }
        }
        System.out.println();
    }

    public static String getReverse(String tag) {
        return new StringBuilder(tag).reverse().toString();
    }
}