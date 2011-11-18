package ru.spbau.bioinf.tagfinder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class UnmatchedStatistics {

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        int[] stat = new int[100];
        int max = 0;
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        int total = 0;
        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            if (!msAlignResults.containsKey(scanId)) {
                List<Peak> peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan));
                GraphUtil.generateEdges(conf, peaks);
                Peak[] bestTag = GraphUtil.findBestTag(peaks);
                int v = bestTag.length - 1;
                //if (v >= 3) {
                System.out.println(scanId + " " + v);
                //}
                stat[v]++;
                if (max < v) {
                    max = v;
                }
                total++;
            }
        }

        System.out.println("\\begin{table}[h]\n" +
                "\\vspace{3mm}\\\n" +
                "{\\centering\n" +
                "\\begin{center}\n" +
                "\\begin{tabular}{|c|c|");

        for (int i = 0; i <= max; i++) {
            System.out.print("c|");
        }

        System.out.println("}\n" +
                "  \\hline\n" +
                "  \\multicolumn{2}{|c|}{} & \\multicolumn{ " + (max + 1) + "}{|c|}{maximum tag length} \\\\\n" +
                "  \\cline{3-" + (max  + 3)+ "}\n" +
                "  \\multicolumn{2}{|c|}{} ");
        for (int i = 0; i <= max; i++) {
            System.out.print(" & " + i);
        }

        System.out.println("\\\\\n" +
                "  \\hline\n" +
                "  \\multirow{2}{*}{spectra} & \\#");

        for (int i = 0; i <= max; i++) {
            System.out.print(" & " + stat[i]);
        }

        System.out.println("\\\\\n" +
                "   & \\%");


        for (int i = 0; i <= max; i++) {
            System.out.print(" & " + ValidTags.df.format(stat[i] * 100d/total));
        }


        System.out.println("\\\\\n" +
                "  \\hline\n" +
                "\\end{tabular}\n" +
                "\\end{center}\n" +
                "\\par}\n" +
                "\\centering\n" +
                "\\caption{The number and percentage of unidentified spectra with a given maximum tag length, for all the observed tag lengths.}\n" +
                "\\vspace{3mm}\n" +
                "\\label{table:unident-tags}\n" +
                "\\end{table}");
    }
}
