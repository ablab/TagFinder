package ru.spbau.bioinf.tagfinder;

import java.io.BufferedReader;
import java.io.File;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class IntencityTableGenerator {

    public static void main(String[] args) throws Exception {
        String file1 = "bar_exp_annotated_proper_none_intencity";
        String file2 = "bar_exp_annotated_correct_none_intencity";

        double[][] res = new double[100][3];
        for (int gap = 1; gap <= 3; gap++) {
            printTable(res, file1, "The number of $(k-d,d)$-spectra in the $ST$ data set for the case of " + gap + "-aa tags, for all $0\\le d\\le k\\le 25$.", gap);
        }
        TexTableGenerator.createThreeRowsTable(res, "Average percentage of proper top-scoring tags", "k", "");

        res = new double[100][3];
        for (int gap = 1; gap <= 3; gap++) {
            printTable(res, file2, "The number of $(k-d,d)$-spectra in the $\\STbar$ data set for the case of " + gap + "-aa tags, for all $0\\le d\\le k\\le 25$.", gap);
        }
        TexTableGenerator.createThreeRowsTable(res, "Average percentage of correct top-scoring tags", "k", "");
    }

    private static void printTable(double[][] res, String file, String caption, int gap) throws Exception {
        BufferedReader in = ReaderUtil.createInputReader(new File("res", "share_" + file + "_" + gap + ".txt"));
        double[] good = new double[100];
        double[] bad = new double[100];
        double[] both = new double[100];

        int max = 0;

        do {
            String s = in.readLine();
            if (s.contains(",")) {
                break;
            }

            if (s.indexOf(" ") == 0) {
                break;
            }
            long[] data = CalculateRelation.getData(s);
            int pos = 2;
            int d = 1;

            while (pos < data.length) {
                max = Math.max(max, d);
                if (data[pos] == 0) {
                    bad[d]++;
                } else if (data[pos + 1] == 0) {
                    good[d]++;
                } else {
                    both[d]++;
                }

                d++;
                pos += 2;
            }
        } while (true);
        double[] ans = new double[max];
        for (int i = 1; i<= max; i++) {
            double total = good[i] + bad[i] + both[i];
            ans[i - 1] = (good[i] + bad[i])*100/total;
        }
        res[gap - 1] = ans;
    }
}
