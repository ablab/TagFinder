package ru.spbau.bioinf.tagfinder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class KDStatistics {

    private Configuration conf;

    private static int gap = 3;

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();
        Map<Integer,Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        KDStatistics kdStat = new KDStatistics(conf);
        Map<KD, Integer> stat = new HashMap<KD, Integer>();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Set<Integer> usedProteins = new HashSet<Integer>();
        Collections.sort(keys);
        System.out.println("gap = " + gap);
        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            if (msAlignResults.containsKey(scanId)) {
                Integer proteinId = msAlignResults.get(scanId);
                if (usedProteins.contains(proteinId)) {
                    //continue;
                }
                usedProteins.add(proteinId);
                if (scanId != 1946) {
                    //continue;
                }
                String sequence = proteins.get(proteinId).getSimplifiedAcids();
                System.out.print(scanId + " ");
                KD kd =
                        kdStat.
                                findKd
                                //findKdBetweenBYAll
                                //findKdBetweenBY
                        (scan, sequence);
                System.out.println(kd.toString() + " " + proteinId);
                if (stat.containsKey(kd)) {
                    stat.put(kd, 1 + stat.get(kd));
                } else {
                    stat.put(kd, 1);
                }
            }
        }
        List<KD> values  = new ArrayList<KD>();
        values.addAll(stat.keySet());
        Collections.sort(values);
        for (KD value : values) {
            System.out.println(value.toString() + " - " + stat.get(value));
        }
    }

    public KDStatistics(Configuration conf) {
        this.conf = conf;
    }

    public KD findKd(Scan scan, String sequence) {
        List<Peak> peaks = scan.createStandardSpectrum();
        return findKd(peaks, sequence);
    }

    public KD findKdBetweenBY(Scan scan, String sequence) {
        List<Peak> peaks = new ArrayList<Peak>();
        peaks.addAll(scan.getPeaks());
        double precursorMassShift = PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan);
        peaks.add(new Peak(0, 0, 0));
        double newPrecursorMass = scan.getPrecursorMass() + precursorMassShift;
        peaks.add(new Peak(newPrecursorMass, 0, 0));
        for (Peak peak : scan.getPeaks()) {
             peaks.add(peak.getYPeak(newPrecursorMass));
        }
        return findKd(peaks, sequence);
    }

    public KD findKdBetweenBYAll(Scan scan, String sequence) {
        List<Double> values = PrecursorMassShiftFinder.getAllPossibleShifts(conf, scan);
        values.add(0d);
        values.add(PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan));

        Collections.sort(values);
        KD ans = new KD(0, 0);
        double prev = -1000;
        int tries = 0;
        for (double precursorMassShift : values) {
            if (precursorMassShift - prev < 0.0001) {
                continue;
            }
            prev = precursorMassShift;
            tries++;
            List<Peak> peaks = new ArrayList<Peak>();
            peaks.addAll(scan.getPeaks());
            for (Peak peak : peaks) {
                peak.clearEdges();
            }
            peaks.add(new Peak(0, 0, 0));
            double newPrecursorMass = scan.getPrecursorMass() + precursorMassShift;
            peaks.add(new Peak(newPrecursorMass, 0, 0));
            for (Peak peak : scan.getPeaks()) {
                 peaks.add(peak.getYPeak(newPrecursorMass));
            }
            KD newKD = findKd(peaks, sequence);
            if (newKD.compareTo(ans) < 0) {
                ans = newKD;
            }
        }
        //System.out.println("tries = " + tries + " " + values.size());
        return ans;
    }

    private KD findKd(List<Peak> peaks, String sequence) {
        Collections.sort(peaks);

        int n = peaks.size();
        for (int i = 0; i < n; i++) {
            peaks.get(i).setComponentId(i);
        }

        GraphUtil.generateGapEdges(conf, peaks, gap);

        boolean done;

        do {
            done = true;
            for (Peak peak : peaks) {
                if (peak.updateComponentId()) {
                    done = false;
                }
            }
        } while (!done);

        int[] kValues = new int[n];
        for (Peak peak : peaks) {
            peak.setMaxPrefix(-1);
        }

        for (Peak peak : peaks) {
            searchK(kValues, 0, peak, new Peak[500] );
        }

        int k = 0;

        for (int kValue : kValues) {
            if (kValue > k) {
                k = kValue;
            }
        }

        //System.out.println("k = " + k);
        kGlobal = k;
        dGlobal = 0;

        String sequenceReversed = new StringBuffer(sequence).reverse().toString();

        for (int i = 0; i < n; i++) {
            Peak peak =  peaks.get(i);
            if (kValues[peak.getComponentId()] == k) {
                getD(peak, sequence);
                getD(peak, sequenceReversed);
            }
        }

        return new KD(k, dGlobal);
    }

    private void searchK(int[] kValues, int len, Peak peak, Peak[] prefix) {
        if (peak.getMaxPrefix() >= len) {
            return;
        }

        prefix[len] = peak;

        int componentId = peak.getComponentId();
        if (kValues[componentId] < len) {
            kValues[componentId] = len;

            //printPrefix(prefix, len);

        }

        for (Peak next : peak.getNext()) {
            searchK(kValues, len + 1, next, prefix);
        }

        peak.setMaxPrefix(len);
    }

    private void printPrefix(Peak[] prefix, int len) {
        if (len < 25) {
            return;
        }
        System.out.print(len + " ");
        for (int i = 0; i < len; i++) {
            Peak p = prefix[i];
            System.out.print(p.getValue() + " ");
            System.out.print(Acid.getAcid(prefix[i + 1].diff(p)).name() + " ");
        }
        System.out.println(prefix[len].getValue());
    }



    private void getD(Peak peak, String sequence) {
        Peak[] prefix = new Peak[500];
        for (int i =0; i < sequence.length() - 1; i++) {
            getD(peak, sequence.substring(i), 0, prefix);
        }
    }

    private void getD(Peak peak, String sequence, int matched, Peak[] prefix) {
        if (dGlobal == kGlobal) {
            return;
        }
        if (matched > dGlobal) {
            dGlobal = matched;
        }

        prefix[matched] = peak;
        int ans = matched;
        for (Peak next : peak.getNext()) {
            double[] limits = conf.getEdgeLimits(peak, next);
            checkNext(sequence, matched, prefix, ans, next, limits);
        }
    }

    private int dGlobal = 0;
    private int kGlobal = 0;

    private int checkNext(String sequence, int matched, Peak[] prefix, int ans, Peak next, double[] limits) {
        double mass = 0;
        for (int i = 0; i < gap; i++) {
            if (sequence.length() <= i) {
                break;
            }
            mass += Acid.getAcid(sequence.charAt(i)).getMass();
            if (limits[0] < mass && limits[1] > mass) {
                getD(next, sequence.substring(i + 1), matched + 1, prefix);
            }
        }
        return ans;
    }
}
