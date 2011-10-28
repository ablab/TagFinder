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

    public static final double TEN_PPM = 0.00001d;

    public static final String VIRTUAL = "virtual";
    public static final String BASIC = "basic";
    public static final String COMBINED = "combined";

    public static final int FULL = 1;
    public static final int BAR = 0;

    private Configuration conf;

    private int gap;
    private static Map<Integer,List<Peak>> msAlignPeaks;
    private Map<KD,Integer> kdStat;
    private int k;
    private int d;
    private Map<Integer,double[]> annotatedSpectrums;

    public ValidTags(Configuration conf) throws Exception {
        this.conf = conf;
        msAlignPeaks = conf.getMSAlignPeaks(conf.getScans());
        annotatedSpectrums = conf.getAnnotatedSpectrums();
    }

    private static NumberFormat df = NumberFormat.getInstance();
    static {
        df.setMaximumFractionDigits(2);
    }

    private boolean addOnes = false;
    private boolean needCorrect = false;

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        ValidTags validTags = new ValidTags(conf);

        for (int gap = 1; gap < 4; gap++) {
            validTags.process(BASIC, FULL, gap, false, false);//I

            validTags.process(BASIC, BAR, gap, false, false);//1

            validTags.process(VIRTUAL, BAR, gap, false, false);//2
            validTags.process(COMBINED, BAR, gap, false, false);
            validTags.process(BASIC, BAR, gap, false, true);

            validTags.process(BASIC, BAR, gap, true, false);//3
            validTags.process(VIRTUAL, BAR, gap, true, false);
        }
    }

    private PrintWriter output;
    private PrintWriter outputKD;


    private void process(String type, int datasetType,  int gap, boolean needCorrect, boolean addOnes) throws Exception {
        this.gap = gap;
        this.needCorrect = needCorrect;
        this.addOnes = addOnes;
        List<Protein> proteins = conf.getProteins();
        Map<Integer,Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        Set<Integer> usedProteins = new HashSet<Integer>();

        double[] global = new double[100];
        int[] count = new int[100];

        String fileName = datasetType == FULL ? "full" : "bar";

        fileName += "_" + type + "_" + gap + "_";
        fileName += needCorrect ? "correct" : "proper";
        if (addOnes) {
            fileName += "_add";
        }
        output = ReaderUtil.createOutputFile(new File("res", "share_" + fileName + ".txt"));
        outputKD = ReaderUtil.createOutputFile(new File("res", "kd_" + fileName + ".txt"));

        System.out.println("fileName = " + fileName);
        kdStat = new HashMap<KD, Integer>();

        //keys.clear(); keys.add(679);
        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();

            if (msAlignResults.containsKey(scanId)) {
                Integer proteinId = msAlignResults.get(scanId);
                if (usedProteins.contains(proteinId) && datasetType == BAR) {
                    continue;
                }
                Protein protein = proteins.get(proteinId);
                long[][] stat = process(scan, protein, type, gap);


                for (int i = 1; i < stat.length; i++) {
                    long good = stat[i][0];
                    long total = good + stat[i][1];
                    if (total > 0) {
                        global[i] += (100d * good) / total;
                        count[i]++;
                    } else {
                        break;
                    }
                }

                printStat(stat, proteinId, scanId);
                usedProteins.add(proteinId);
            }
        }
        for (int i = 1; i < count.length; i++) {
            int n = count[i];
            if (n > 0) {
                String text = df.format(global[i] / n) + " ";
                System.out.print(text);
                output.print(text);
            } else {
                break;
            }
        }
        List<KD> values  = new ArrayList<KD>();
        values.addAll(kdStat.keySet());
        Collections.sort(values);
        for (KD value : values) {
            String text = value.toString() + " - " + kdStat.get(value);
            outputKD.println(text);
        }

        output.close();
        outputKD.close();
        System.out.println();
    }

    public long[][] process(Scan scan, Protein protein, String type, int gap) {
        wrongCache.clear();
        correctCache.clear();
        this.gap = gap;
        String sequence = protein.getSimplifiedAcids();
        List<Peak> peaks;

        int scanId = scan.getId();
        if (VIRTUAL.equals(type)) {
            peaks = msAlignPeaks.get(scanId);
        } else {
            PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan);
            peaks = scan.createSpectrumWithYPeaks(0);
        }
        if (addOnes) {
            peaks = addOnes(peaks);
        }

        GraphUtil.generateGapEdges(conf, peaks, gap);
        if (BASIC.equals(type)) {
            filterMonotags(peaks);
        }

        int proteinId = protein.getProteinId();
        double[] proteinSpectrum = needCorrect ? annotatedSpectrums.get(scanId): ShiftEngine.getSpectrum(sequence);
        SpectrumResult spectrumResult = printGappedTagInfo(peaks, proteinSpectrum, scan.getPrecursorMass());
        spectrumResult.output(kdStat, outputKD, scanId, proteinId);

        return spectrumResult.stat;
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
                    if (peak.getIntensity() != 0 && next.getIntensity() != 0) {
                        iterator.remove();
                    }
                }
            }
        }
    }

    private SpectrumResult printGappedTagInfo(List<Peak> peaks, double[] proteinSpectrum, double precursorMass) {
        long[][] stat = new long[1000][2];
        KD kd = new KD(0, 0);
        List<List<Peak>> components = GraphUtil.getComponentsFromGraph(peaks);
        for (List<Peak> component : components) {
            k = 0;
            d = 0;
            for (Peak peak : component) {
                Set<Integer> starts = new HashSet<Integer>();
                for (int i = 0; i < proteinSpectrum.length; i++) {
                    if (needCorrect) {
                        double mass = proteinSpectrum[i];
                        double limit = peak.getPeakType() == PeakType.Y ? TEN_PPM * 1.5 * precursorMass : TEN_PPM * peak.getMass();
                        if (Math.abs(mass - peak.getValue()) < limit) {
                            starts.add(i);
                        }
                    } else {
                        starts.add(i);
                    }
                }
                processGappedTags(stat, peak, 0, proteinSpectrum, starts);
            }
            KD newKD = new KD(k, d);
            if (newKD.compareTo(kd) < 0) {
                kd = newKD;
            }
        }

        return new SpectrumResult(stat, kd);
    }

    private Map<Peak, List<CorrectResult>> correctCache = new HashMap<Peak, List<CorrectResult>>();

    private void processGappedTags(long[][] stat, Peak peak, int prefix, double[] proteinSpectrum, Set<Integer> starts) {
        List<CorrectResult> results = correctCache.get(peak);
        if (results != null) {
            for (CorrectResult result : results) {
                if (result.check(starts)) {
                    long[] deltaGood = result.getDeltaGood();
                    long[] deltaBad = result.getDeltaBad();
                    for (int i = 0; i < deltaGood.length; i++) {
                        int pos = prefix + i;
                        stat[pos][0] += deltaGood[i];
                        stat[pos][1] += deltaBad[i];
                        if (deltaGood[i] > 0) {
                            if (pos > k) {
                                k = pos;
                            }
                            if (pos > d) {
                                d = pos;
                            }
                        }
                        if (deltaBad[i] > 0) {
                            if (pos > k) {
                                k = pos;
                            }
                        }
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
            if (prefix > k) {
                k = prefix;
            }
            if (prefix > d) {
                d = prefix;
            }
        }

        List<Peak> nextPeaks = peak.getNext();
        for (Peak next : nextPeaks) {
            double[] limits = conf.getEdgeLimits(peak, next);
            Set<Integer> nextStarts = getNextStarts(proteinSpectrum, starts, limits, gap);
            if (nextStarts.size() == 0) {
                processWrongGappedTags(stat, next, prefix + 1);
            } else {
                processGappedTags(stat, next, prefix + 1, proteinSpectrum, nextStarts);
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
        correctCache.get(peak).add(new CorrectResult(starts, deltaGood, deltaBad));
    }

    private Map<Peak, long[]> wrongCache = new HashMap<Peak, long[]>();
    private void processWrongGappedTags(long[][] stat, Peak peak, int prefix) {
        if (wrongCache.containsKey(peak)) {
            long[] delta = wrongCache.get(peak);
            for (int i = 0; i < delta.length; i++) {
                int pos = prefix + i;
                stat[pos][1] += delta[i];
                if (delta[i] > 0) {
                    if (pos > k) {
                        k = pos;
                    }
                }
            }
            return;
        }

        long[] old = new long[50];
        for (int i = 0; i < old.length; i++) {
            old[i] = stat[i + prefix][1];
        }
        List<Peak> nextPeaks = peak.getNext();
        stat[prefix][1]++;
        if (prefix > k) {
            k = prefix;
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

    public static Set<Integer> getNextStarts(double[] proteinSpectrum, Set<Integer> starts, double[] limits, int gap) {
        Set<Integer> nextStarts = new HashSet<Integer>();
        for (int pos : starts) {
            for (int i = 1; i <= gap; i++) {
                if (pos + i < proteinSpectrum.length ) {
                    double m = proteinSpectrum[pos + i] - proteinSpectrum[pos];
                    if (limits[0] < m && m < limits[1]) {
                        nextStarts.add(pos + i);
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
        private long[] deltaGood;
        private long[] deltaBad;

        public CorrectResult(Set<Integer> starts, long[] deltaGood, long[] deltaBad) {
            this.starts = starts;
            this.deltaGood = deltaGood;
            this.deltaBad = deltaBad;
        }

        public boolean check(Set<Integer> startsNew) {
            return starts.equals(startsNew);
        }

        public long[] getDeltaGood() {
            return deltaGood;
        }

        public long[] getDeltaBad() {
            return deltaBad;
        }
    }

    public static class SpectrumResult {
        private long[][] stat;
        private KD kd;

        public SpectrumResult(long[][] stat, KD kd) {
            this.stat = stat;
            this.kd = kd;
        }

        public void output(Map<KD,Integer> kdStat, PrintWriter outputKD, int scanId, int proteinId) {
            if (kdStat.containsKey(kd)) {
                kdStat.put(kd, 1 + kdStat.get(kd));
            } else {
                kdStat.put(kd, 1);
            }
            outputKD.println(scanId + " " + kd.toString() + " - " + proteinId);
        }
    }
}