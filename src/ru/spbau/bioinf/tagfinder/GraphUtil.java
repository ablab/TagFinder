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

    public static Peak[] findBestTag(List<Peak> peaks) {
        initMaxPrefix(peaks);
        Peak[] bestTag;
        bestTag = new Peak[]{};
        for (Peak peak : peaks) {
                bestTag = findBestTag(peak, bestTag, 0, new Peak[500]);
            }
        return bestTag;
    }

    public static void initMaxPrefix(List<Peak> peaks) {
        for (Peak peak : peaks) {
            peak.setMaxPrefix(-1);
        }
    }

    public static Peak[] findBestTag(Peak peak, Peak[] best, int len, Peak[] prefix) {
        if (peak.getMaxPrefix() >= len) {
            return best;
        }

        prefix[len] = peak;

        if (len >= best.length) {
            best = new Peak[len +1];
            System.arraycopy(prefix, 0, best, 0, len + 1);
        }

        for (Peak next : peak.getNext()) {
            best = findBestTag(next, best, len + 1, prefix);
        }

        peak.setMaxPrefix(len);

        return best;
    }
}
