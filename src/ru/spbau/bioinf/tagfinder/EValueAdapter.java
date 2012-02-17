package ru.spbau.bioinf.tagfinder;

import edu.ucsd.msalign.align.MsAlign;
import edu.ucsd.msalign.align.PropertyUtil;
import edu.ucsd.msalign.align.idevalue.IdEValue;
import edu.ucsd.msalign.align.prsm.PrSM;
import edu.ucsd.msalign.seq.MsAlignSeq;
import java.util.List;
import java.util.Map;
import java.util.Properties;

public class EValueAdapter {

    private static MsAlign align;

    public static String getVersion() {
        return "zero0.7";
    }

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        init(conf);
        Map<Integer, Scan> scans = conf.getScans();
        conf.getMSAlignResults();
        for (int i = 0; i < input.length; i += 2) {
            Scan scan = scans.get(input[i]);
            int proteinId = input[i + 1];
            System.out.println("Processing scan " + scan.getName() + " protein " + proteinId);
            output(calculateEValue(scan, proteinId));
        }
    }

    public static void init(Configuration conf) throws Exception {
        conf.getScans();
        Properties prop = PropertyUtil.getDefaultProperties();
        prop.put("databaseFileName", conf.getProteinDatabaseFile().getAbsolutePath());
        prop.put("spectrumFileName", conf.getMsalignFile().getAbsolutePath());
        prop.put("activation", "CID");
        prop.put("cysteineProtection", "C57");
        align = new MsAlign(prop);
        //align.prepare();
        align.readSpectra();

    }


    public static synchronized List<PrSM[]> calculateEValue(Scan scan, int... proteinIds) throws Exception {
        List<PrSM[]> prsms;
        synchronized (align) {
            align.updateSeqs(proteinIds);
            prsms = align.getPrsms(scan.getId());
        }
        return prsms;
    }

    private static void output(List<PrSM[]> prsms) {
        if (prsms != null) {
            for (PrSM[] prsmA : prsms) {
                for (PrSM prsm : prsmA) {
                    System.out.println("Shift " + prsm.getInterShiftNum() + " Alignment type " + prsm.getEValueAlignType() + " score " + prsm.getUniqueScr() + " evalue " + prsm.getEValue());
                }
            }
        }
    }

    public static PrSM getBestEValue(Scan scan, int proteinId) throws Exception {
        List<PrSM[]> prsms = calculateEValue(scan, proteinId);
        double best = 9E100;
        PrSM ans = null;
        if (prsms != null) {
            for (PrSM[] prsmA : prsms) {
                for (PrSM prsm : prsmA) {
                    double newValue = prsm.getEValue();
                    if (newValue < best) {
                        best = newValue;
                        ans = prsm;
                    }
                }
            }
        }
        return ans;
    }

    private static int[] input = new int[]{
            1082, 4375,
            1082, 392,
            1082, 1007,
            1082, 1573,
            1082, 2508,
            1082, 3392,
            1082, 3874,
            1082, 3949,

            1083, 4149,
            1083, 3949,

            1242, 1100,
            1242, 3302,
            1242, 3393,

            1251, 1100,
            1251, 1618,
            1251, 2013,
            1251, 3302,

            1252, 3302,

            1274, 3302,
            1274, 471,
            1274, 4098,

            1322, 3091,
            1322, 1453,

            1323, 2002,
            1323, 1453,

            1326, 152,
            1326, 1453,

            1375, 3307,
            1375, 381,
            1375, 1277,

            1378, 2272,
            1378, 2115,
            1378, 3734,
            1378, 4098,

            1711, 3311,
            2167, 4368,
            2167, 492
    };
}
