package ru.spbau.bioinf.tagfinder;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class KDStatistics {

    private static Configuration conf;

    public static void main(String[] args) throws Exception {
        conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();
        Map<Integer,Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        KDStatistics kdStat = new KDStatistics();
        Map<KD, Integer> stat = new HashMap<KD, Integer>();
        for (Scan scan : scans.values()) {
            int scanId = scan.getId();
            if (msAlignResults.containsKey(scanId)) {
                String sequence = proteins.get(msAlignResults.get(scanId)).getSimplifiedAcids();
                KD kd = kdStat.findKd(scan, sequence);
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

    public KD findKd(Scan scan, String sequence) {
        KD ans = null;

        List<Peak> peaks = scan.getPeaks();
        Collections.sort(peaks);
        double[] edges = new double[Acids.acids.size()];
        int cur = 0;
        for (double v : Acids.acids.values()) {
            edges[cur] = v;
            cur++;
        }

        int n = peaks.size();
        for (int i = 0; i < n; i++) {
            peaks.get(i).setComponentId(i);
        }

        for (int i = 0; i < n; i++) {
            Peak peak = peaks.get(i);
            for (int j = i+1; j < n; j++) {
                Peak next =  peaks.get(j);
                double diff = next.diff(peak);
                double error =  next.getValue() * conf.getPpmCoef();
                double min = diff - error;
                double max = diff + error;
                for (double edge : edges) {
                    if (min < edge && max > edge) {
                        peak.addNext(next);
                        break;
                    }
                }
            }
        }

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
            peak.setProcessed(false);
        }

        for (Peak peak : peaks) {
            searchK(kValues, 0, peak);
        }

        int maxK = 0;

        for (int kValue : kValues) {
            if (kValue > maxK) {
                maxK = kValue;
            }
        }

        ans = new KD(maxK, 0);

        return ans;
    }

    private void searchK(int[] kValues, int len, Peak peak) {
        if (peak.isProcessed()) {
            return;
        }

        int componentId = peak.getComponentId();
        if (kValues[componentId] <= len) {
            kValues[componentId] = len;
        }

        for (Peak next : peak.getNext()) {
            searchK(kValues, len + 1, next);
        }

        peak.setProcessed(true);
    }
}
