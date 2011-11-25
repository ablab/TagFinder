package ru.spbau.bioinf.tagfinder;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class YFinder {

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();
        Map<Integer,Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        //keys.clear();keys.add(3008);
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        int good = 0;
        for (int scanId : keys) {
            if (scanId < 502) {
                continue;
            }
            List<Protein> bestProteins = new ArrayList<Protein>();
            int bestScore = 1;
            double bestShift = 0;
            List<Peak> peaks = scans.get(scanId).getPeaks();
            Collections.sort(peaks);
            double[] pw = new double[peaks.size()];
            for (int i = 0; i < pw.length; i++) {
                pw[i] = peaks.get(i).getMass()
                        - Consts.WATER
                ;
            }

            double[] pk = new double[peaks.size()];
            for (int i = 0; i < pw.length; i++) {
                pk[i] = peaks.get(i).getMass()
                        + Acid.K.getMass() - 30.0486;
            }

            double[] p = new double[peaks.size()];
            for (int i = 0; i < p.length; i++) {
                p[i] = peaks.get(i).getMass()
                //        - Consts.WATER
                ;
            }

            long t = 0;

            for (Protein protein : proteins) {
                double[] yEnds = protein.getYEnds();
                double[] bEnds = protein.getBEnds();
                if (yEnds.length < 1 || pw.length  <1) {
                    continue;
                }
                if (protein.getProteinId() != 3941) {
                    //continue;
                }

                /*
                double yscore = 3 * ShiftEngine.getScore(pw, yEnds, 0) + ShiftEngine.getScore(p, yEnds, 0);
                double kscore = 3 * ShiftEngine.getScore(pk, yEnds, 0);
                int score = (int)Math.round(Math.max(yscore, kscore)
                        + 3 * ShiftEngine.getScore(p, bEnds, 0) + ShiftEngine.getScore(pw, bEnds, 0)
                );
                */
                if (protein.getProteinId() == 570) {
                    //System.out.println();
                } else {
                    //continue;
                }
                long start = System.currentTimeMillis();
                double shift = ShiftEngine.getBestShiftFast(p, yEnds);
                long bst = System.currentTimeMillis();
                t += (bst - start);
                int score = (int)ShiftEngine.getScore(p, yEnds, shift);

                if (score > bestScore) {
                    bestProteins.clear();
                    bestScore = score;
                    bestShift = shift;
                }
                if (score >= bestScore) {
                    bestProteins.add(protein);
                    //System.out.println(scanId + " " + protein.getProteinId() + " "  + score);
                }
            }
            //if (bestScore > 3) {
                Integer proteinId = msAlignResults.get(scanId);
                System.out.print(t + " " + scanId + " " + proteinId + " " + bestScore + " " + bestShift + " " );

                if (proteinId != null) {
                    String res = "fail";
                    for (Protein protein : bestProteins) {
                        if (protein.getProteinId() == proteinId.intValue()) {
                            res = "good";
                        }
                    }
                    System.out.print(res + " ");
                    if (res.equals("fail")) {
                        //System.out.println(scanId + " " + proteinId + " " + ValidTags.getReverse(proteins.get(proteinId).getSimplifiedAcids()));
                        //System.out.print(scanId + " " + proteinId + " ");
                        if (bestProteins.size() == 1) {
                            for (Protein protein : bestProteins) {
                                //System.out.print(protein.getProteinId() + " ");
                            }
                        }
                    } else {
                        good++;
                    }
                    System.out.print(" " + conf.getEvalues().get(scanId) + " ");
                }


            if (bestScore > 1) {
                for (Protein protein : bestProteins) {
                    System.out.print(protein.getProteinId() + " ");
                }
            }
            System.out.println();
                //System.out.println(bestScore);
            //}
            //System.out.println(scanId);
        }
        System.out.println("good = " + good);
    }
}
