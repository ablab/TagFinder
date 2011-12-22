package ru.spbau.bioinf.tagfinder;

import java.util.Collections;
import java.util.List;

public class PeaksFilter {

    public static List<Peak> getPeaks(Configuration conf, List<Peak> peaks) {
        Collections.sort(peaks);
        for (int i = 0; i < peaks.size() - 1; i++) {
            Peak p1 = peaks.get(i);
            Peak p2 = peaks.get(i + 1);
            if (Math.abs(p1.getValue() - p2.getValue()) < 0.1) {
                boolean isSame = true;
                for (Peak p3 : peaks) {
                    for (Acid acid : Acid.values()) {
                        if (acid.match(conf.getEdgeLimits(p2, p3)) && !acid.match(conf.getEdgeLimits(p1, p3))) {
                            isSame = false;
                            break;
                        }
                    }
                    if (!isSame) {
                        break;
                    }
                }
                if (!isSame) {
                    continue;
                }

                peaks.remove(p2);
                peaks.remove(p1);
                double mass = (p1.getMass() + p2.getMass()) / 2;
                peaks.add(new Peak(mass, p1.getIntensity(), p1.getCharge()));
                return getPeaks(conf, peaks);
            }
        }
        return peaks;
    }
}
