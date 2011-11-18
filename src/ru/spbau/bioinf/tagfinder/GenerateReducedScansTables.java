package ru.spbau.bioinf.tagfinder;

import java.io.BufferedReader;
import java.io.File;
import java.util.HashMap;
import java.util.Map;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class GenerateReducedScansTables {

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args, UnmatchedScansGenerator.SHARED_MODE);
        ValidTags validTags = new ValidTags(conf);
        ValidTags.parentDir = "res-reduced";
        validTags.process(ValidTags.INPUT_EXP, ValidTags.TARGET_BASE, ValidTags.MATCH_PROPER, ValidTags.BY_NONE, ValidTags.FULL, 1, false, false);
        Map<Integer, String> oldKd = new HashMap<Integer, String>();
        BufferedReader inOld = ReaderUtil.createInputReader(new File("res", "kd_full_exp_base_proper_none_1.txt"));
        do {
            String s = inOld.readLine();
            if (s == null) {
                break;
            }
            if (s.startsWith("(")) {
                continue;
            }
            String[] data = s.split("\\D+");
            int scanId = Integer.parseInt(data[0]);
            int k = Integer.parseInt(data[1]);
            int d = Integer.parseInt(data[2]);
            oldKd.put(scanId, "(" + k + "," + d + ")");
        } while (true);

        int total = 0;
        int[] stat = new int[100];
        int maxK = 0;

        BufferedReader in = ReaderUtil.createInputReader(new File("res-reduced", "kd_full_exp_base_proper_none_1.txt"));

        System.out.println("\\begin{table}[h]\\footnotesize\n" +
                "\\vspace{3mm}\\\n" +
                "{\\centering\n" +
                "\\begin{center}\n" +
                "\\begin{tabular}{|c|c|c|c|}\n" +
                "  \\hline\n" +
                "  \\multicolumn{2}{|c|}{PrSM by MS-Align+} & \\multirow{2}{*}{$(k,d)$} & maximum tag length \\\\\n" +
                "  \\cline{1-2}\n" +
                "  scan & protein & &  in the reduced spectrum \\\\\n" +
                "  \\hline");
        do {
            String s = in.readLine();
            if (s == null) {
                break;
            }
            if (s.startsWith("(")) {
                continue;
            }
            String[] data = s.split("\\D+");
            int scanId = Integer.parseInt(data[0]);
            int k = Integer.parseInt(data[1]);
            int d = Integer.parseInt(data[2]);
            int proteinId = Integer.parseInt(data[3]);
            if (k >= 3) {
                System.out.println(scanId + "&  " + proteinId + " & " + oldKd.get(scanId) + " & " + k + "  \\\\");
                total++;
                stat[k]++;
                if (k > maxK) {
                    maxK = k;
                }
            }
        } while (true);

        System.out.println("  \\hline\n" +
                "\\end{tabular}\n" +
                "\\end{center}\n" +
                "\\par}\n" +
                "\\centering\n" +
                "\\caption{Candidate mixed spectra.}\n" +
                "\\vspace{3mm}\n" +
                "\\label{table:mixed-cand}\n" +
                "\\end{table}");

        System.out.println("\\begin{table}[h]\n" +
                "\\vspace{3mm}\\\n" +
                "{\\centering\n" +
                "\\begin{center}\n" +
                "\\begin{tabular}{|");
        for (int i = 1; i <= maxK; i++) {
            System.out.print("c|");
        }

        System.out.println("}\n" +
                "  \\hline\n" +
                "  \\multicolumn{2}{|c|}{} & \\multicolumn{" + (maxK - 3 + 1)+ "}{|c|}{maximum tag length} \\\\\n" +
                "  \\cline{3- " + maxK + "}\n" +
                "  \\multicolumn{2}{|c|}{} ");
        for (int i = 3; i <= maxK; i++) {
            System.out.print(" & " + i);
        }

        System.out.println(" \\\\\n" +
                "  \\hline\n" +
                "  \\multirow{2}{*}{spectra} & \\# ");
        for (int i = 3; i <= maxK; i++) {
            System.out.print(" & " + stat[i]);
        }
        System.out.println(" \\\\\n" +
                "   & \\%");

        for (int i = 3; i <= maxK; i++) {
            System.out.print(" & " + ValidTags.df.format(stat[i]*100d/total));
        }

        System.out.println("\\\\\n" +
                "  \\hline\n" +
                "\\end{tabular}\n" +
                "\\end{center}\n" +
                "\\par}\n" +
                "\\centering\n" +
                "\\caption{The number and percentage of candidate mixed spectra with the maximum tag length in the reduced spectra at least $3$.}\n" +
                "\\vspace{3mm}\n" +
                "\\label{table:mixed-stat}\n" +
                "\\end{table}");

    }
}
