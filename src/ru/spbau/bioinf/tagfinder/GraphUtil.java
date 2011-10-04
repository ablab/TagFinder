package ru.spbau.bioinf.tagfinder;

import java.util.List;

public class GraphUtil {
    public static void generateEdges(Configuration conf, List<Peak> peaks) {
        int n = peaks.size();
        for (int i = 0; i < n; i++) {
            Peak peak = peaks.get(i);
            for (int j = i+1; j < n; j++) {
                Peak next =  peaks.get(j);
                double[] limits = conf.getEdgeLimits(peak, next);
                for (Acid acid : Acid.values()) {
                    if (acid.match(limits)) {
                        peak.addNext(next);
                        break;
                    }
                }
            }
        }
    }
}
