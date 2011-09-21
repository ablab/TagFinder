package ru.spbau.bioinf.tagfinder.util;

import java.util.Collections;
import java.util.List;

public class StatisticsUtil {

    public static double getAverage(List<Double> values) {
        Collections.sort(values);
        int n = values.size();
        if (n == 0) {
            return 0;
        }
        int cut = n /10;
        double total = 0;
        for (int i = cut; i < values.size() - cut; i++) {
            total += values.get(i);
        }
        return total / (n - 2 * cut);
    }

}
