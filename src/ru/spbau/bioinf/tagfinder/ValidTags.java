package ru.spbau.bioinf.tagfinder;

import java.io.File;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class ValidTags {

    public static final String VIRTUAL = "virtual";
    public static final String STANDARD = "standard";
    public static final String REFLEX = "reflex";

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

    private boolean needReverseTag = false;
    private boolean monoTagsOnly = false;
    private boolean addOnes = false;
    private int maxTagLength = 90;

    public void setMonoTagsOnly(boolean monoTagsOnly) {
        this.monoTagsOnly = monoTagsOnly;
    }

    public void setAddOnes(boolean addOnes) {
        this.addOnes = addOnes;
    }

    public void setNeedReverseTag(boolean needReverseTag) {
        this.needReverseTag = needReverseTag;
    }

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        ValidTags validTags = new ValidTags(conf);
        String type = REFLEX;
        //validTags.setMonoTagsOnly(true);
        //validTags.setAddOnes(true);
        validTags.process(type, 1);
        validTags.process(type, 2);
        validTags.process(type, 3);
    }

    private PrintWriter output;

    private void process(String type, int gap) throws Exception {
        process(type, gap, 100);
    }

    private void process(String type, int gap, int maxTagLength) throws Exception {
        this.gap = gap;
        this.maxTagLength = maxTagLength;
        List<Protein> proteins = conf.getProteins();
        Map<Integer,Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        Set<Integer> usedProteins = new HashSet<Integer>();
        msAlignPeaks = conf.getMSAlignPeaks();

        double[] global = new double[100];
        int[] count = new int[100];


        setNeedReverseTag(STANDARD.equals(type));
        System.out.println("type = " + type);
        System.out.println("gap = " + this.gap);
        System.out.println("monoTagsOnly = " + monoTagsOnly);
        System.out.println("addOnes = " + addOnes);
        System.out.println("needReverseTag = " + needReverseTag);

        String fileName = type + "_" + gap;
        if (monoTagsOnly) {
            fileName += "_mono";
        }
        if (addOnes) {
            fileName += "_add";
        }
        output = ReaderUtil.createOutputFile(new File("vt", fileName + ".txt"));
        //keys.clear();keys.add(502);
        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            if (msAlignResults.containsKey(scanId)) {
                Integer proteinId = msAlignResults.get(scanId);
                if (usedProteins.contains(proteinId)) {
                    continue;
                }
                usedProteins.add(proteinId);
                Protein protein = proteins.get(proteinId);
                long[][] stat = process(scan, protein, type, gap);
                for (int i = 1; i < stat.length; i++) {
                    long good = stat[i][0];
                    long total = good + stat[i][1];
                    if (total > 0) {
                        global[i] += (100d * good)/total;
                        count[i]++;
                    } else {
                        break;
                    }
                }

                printStat(stat, proteinId, scanId);

                //break;
            }
        }
        output.close();
        for (int i = 1; i < count.length; i++) {
            int n = count[i];
            if (n > 0) {
                System.out.print( df.format(global[i]/n) + " ");
            } else {
                break;
            }
        }
        System.out.println();
    }

    public long[][] process(Scan scan, Protein protein, String type, int gap) {
        wrongCache.clear();
        correctCache.clear();
        this.gap = gap;
        String sequence = protein.getSimplifiedAcids();
        List<Peak> peaks = null;

        if (VIRTUAL.equals(type)) {
            peaks = msAlignPeaks.get(scan.getId());
        } else if (STANDARD.equals(type)) {
            peaks = scan.createStandardSpectrum();
        } else if (REFLEX.equals(type)) {
            peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan));
        } else {
            throw new RuntimeException("Wrong spectrum type " +  type);
        }

        if (addOnes) {
            peaks = addOnes(peaks);
        }

        GraphUtil.generateGapEdges(conf, peaks, gap);
        if (monoTagsOnly) {
                filterMonotags(peaks);
        }

        return printGappedTagInfo(peaks, sequence, getReverse(sequence));

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

    private long[][] printGappedTagInfo(List<Peak> peaks, String sequence, String reverseSequence) {
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
        return stat;
    }

    private Map<Peak, List<CorrectResult>> correctCache = new HashMap<Peak, List<CorrectResult>>();

    private void processGappedTags(long[][] stat, Peak peak, int prefix, String sequence, String reverseSequence, Set<Integer> starts, Set<Integer> reverseStarts) {
        List<CorrectResult> results = correctCache.get(peak);
        if (results != null) {
            for (CorrectResult result : results) {
                if (result.check(starts, reverseStarts)) {
                    long[] deltaGood = result.getDeltaGood();
                    long[] deltaBad = result.getDeltaBad();
                    for (int i = 0; i < deltaGood.length; i++) {
                        stat[prefix + i][0] += deltaGood[i];
                        stat[prefix + i][1] += deltaBad[i];
                    }
                    return;
                }
            }
        }
        long[] oldGood = new long[50];
        for (int i = 0; i < oldGood.length; i++) {
            oldGood[i] = stat[i + prefix][0];
        }

        long[] oldBad = new long[50];
        for (int i = 0; i < oldBad.length; i++) {
            oldBad[i] = stat[i + prefix][1];
        }

        if (prefix > 0) {
            stat[prefix][0]++;
        }
        if (prefix > maxTagLength) {
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

        long[] deltaGood = new long[50];
        long[] deltaBad = new long[50];
        for (int i = 0; i < oldGood.length; i++) {
            deltaGood[i] = stat[i + prefix][0] - oldGood[i];
            deltaBad[i] = stat[i + prefix][1] - oldBad[i];
        }
        if (!correctCache.containsKey(peak)) {
            correctCache.put(peak, new ArrayList<CorrectResult>());
        }
        correctCache.get(peak).add(new CorrectResult(starts, reverseStarts, deltaGood, deltaBad));
    }

    private Map<Peak, long[]> wrongCache = new HashMap<Peak, long[]>();
    private void processWrongGappedTags(long[][] stat, Peak peak, int prefix) {
        if (wrongCache.containsKey(peak)) {
            long[] delta = wrongCache.get(peak);
            for (int i = 0; i < delta.length; i++) {
                stat[prefix + i][1] += delta[i];
            }
            return;
        }

        long[] old = new long[50];
        for (int i = 0; i < old.length; i++) {
            old[i] = stat[i + prefix][1];
        }
        List<Peak> nextPeaks = peak.getNext();
        stat[prefix][1]++;

        if (prefix > maxTagLength) {
            return;
        }
        for (Peak next : nextPeaks) {
            processWrongGappedTags(stat, next, prefix + 1);
        }
        long[] delta = new long[50];
        for (int i = 0; i < old.length; i++) {
            delta[i] = stat[i + prefix][1] - old[i];
        }
        wrongCache.put(peak, delta);
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
            if (gap > 1) {
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

        }
        return nextStarts;
    }

    private void printStat(long[][] stat, Integer proteinId, int scanId) {
        output.print(scanId + " " +  proteinId);
        for (int i = 1; i < stat.length; i++) {
            long good = stat[i][0];
            long total = good + stat[i][1];
            if (total > 0) {
                //output.print(" " + df.format((100d * good)/total) + " " + stat[i][0] + " " + stat[i][1] + " ");
                output.print(" " + stat[i][0] + " "  + stat[i][1]);
            } else {
                break;
            }
        }
        output.println();
        output.flush();
    }

    public static String getReverse(String tag) {
        return new StringBuilder(tag).reverse().toString();
    }

    public static class CorrectResult {
        private Set<Integer> starts;
        private Set<Integer> reverseStarts;
        private long[] deltaGood;
        private long[] deltaBad;

        public CorrectResult(Set<Integer> starts, Set<Integer> reverseStarts, long[] deltaGood, long[] deltaBad) {
            this.starts = starts;
            this.reverseStarts = reverseStarts;
            this.deltaGood = deltaGood;
            this.deltaBad = deltaBad;
        }

        public boolean check(Set<Integer> startsNew, Set<Integer> reverseStartsNew ) {
            return starts.equals(startsNew) && reverseStarts.equals(reverseStartsNew);
        }

        public long[] getDeltaGood() {
            return deltaGood;
        }

        public long[] getDeltaBad() {
            return deltaBad;
        }
    }
}