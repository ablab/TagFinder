package ru.spbau.bioinf.palign;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import ru.spbau.bioinf.evalue.EValueServer;
import ru.spbau.bioinf.tagfinder.Acid;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.Consts;
import ru.spbau.bioinf.tagfinder.Peak;
import ru.spbau.bioinf.tagfinder.Protein;
import ru.spbau.bioinf.tagfinder.Scan;

public class Aligner {

    private static List<Protein> proteins;
    private static Map<Integer,Scan> scans;

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        System.out.println("MS-Align results: " + msAlignResults.keySet().size());
        proteins = conf.getProteins();
        scans = conf.getScans();
        EValueServer.init(args);
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        int good = 0;
        int bad = 0;
        for (int scanId : keys) {
            if (scans.get(scanId).getPeaks().size() == 0) {
                continue;
            }
            System.out.print(scanId + " ");
            int bestScore = 0;
            int bestProtein = 0;
            for (Protein protein : proteins) {
                int score = findAlignment(scanId, protein.getProteinId());
                if (score > bestScore) {
                    bestScore = score;
                    bestProtein = protein.getProteinId();
                }
            }
            double evalue = EValueServer.getEvalue(scanId, bestProtein);
            if (evalue < 0.0024) {
                good++;
            } else {
                bad++;
            }

            System.out.println(bestProtein + " " + bestScore + " " + evalue + " " + good + " " + bad);
        }
    }

    public static int findAlignment(int scanId, int proteinId) {
        Scan scan = scans.get(scanId);
        Protein protein = proteins.get(proteinId);
        Alignment alignment = findAlignment(scan, protein);

        return alignment.getCleavages().get(0).getSupports().size();
    }

    public static Alignment findAlignment(Scan scan, Protein protein) {
        Alignment alignment = new Alignment(scan, protein);
        List<Peak> peaks = scan.getPeaks();
        Collections.sort(peaks);
        double[] masses = new double[peaks.size()];
        for (int i = 0; i < masses.length; i++) {
            masses[i] = peaks.get(i).getMass();
        }
        String sequence = protein.getSimplifiedAcids();

        List<Double> diffsB = getDiffs(masses, sequence);
        Cleavage cleavage = getCleavage(peaks, sequence, bestShift(diffsB), false);
        alignment.addCleavage(cleavage);

        String reverseSequence = new StringBuilder(sequence).reverse().toString();
        List<Double> diffsY = getDiffs(masses, reverseSequence);
        cleavage = getCleavage(peaks, reverseSequence, bestShift(diffsY), true);
        alignment.addCleavage(cleavage);

        return alignment;
    }

    private static Cleavage getCleavage(List<Peak> peaks, String sequence, double center, boolean isReverse) {
        double prefix = 0;
        int cleavagePos = 0;
        double proteinModification = center;
        int len = sequence.length();
        for (int i = 0; i < len; i++) {

            Acid acid = Acid.getAcid(sequence.charAt(i));
            if (acid == null) {
                continue;
            }
            prefix += acid.getMass();
            if (Math.abs(center + prefix) < Math.abs(proteinModification)) {
                proteinModification = center + prefix;
                cleavagePos = i + 1;
            }

        }
        if (Math.abs(proteinModification) < 0.1) {
            proteinModification = 0;
        }

        Cleavage cleavage = new Cleavage(isReverse ? len - cleavagePos : cleavagePos, isReverse ? - 1 : 1, proteinModification);
        double peptideMass = 0;
        for (int pos = cleavagePos; pos < len; pos++) {
            peptideMass += Acid.getAcid(sequence.charAt(pos)).getMass();
            PeptideSupport support = null;
            for (Peak peak : peaks) {
                for (ModificationType modificationType : ModificationType.values()) {
                    double error = modificationType.getError(peak, peptideMass + proteinModification);
                    if (Math.abs(error) < 0.2) {
                        if (support == null) {
                            support = new PeptideSupport(isReverse ? len - pos - 1 : pos + 1);
                        }
                        support.addModification(new PeptideModification(peak, modificationType, error));
                    }
                }
            }
            if (support != null) {
                cleavage.addSupport(support);
            }
        }
        return cleavage;
    }

    private static double bestShift(List<Double> diffsB) {
        List<LinkedList<Double>> clusters = new ArrayList<LinkedList<Double>>();
        LinkedList<Double> nonsense = new LinkedList<Double>();
        nonsense.add(-10E100d);
        clusters.add(nonsense);
        for (Double diff : diffsB) {
            LinkedList<Double> last = clusters.get(clusters.size() - 1);
            if (diff - last.getLast() < 0.1) {
                last.add(diff);
            } else {
                LinkedList<Double> cluster = new LinkedList<Double>();
                cluster.add(diff);
                clusters.add(cluster);
            }
        }

        double[] centers = new double[clusters.size()];

        for (int i = 0; i < centers.length; i++) {
            LinkedList<Double> cluster = clusters.get(i);
            double v = 0;
            for (Double diff : cluster) {
                v += diff;
            }
            centers[i] = v / cluster.size();
        }


        int bestScore = 0;
        int bestClusterId = 0;

        for (int i = 1; i < clusters.size(); i++) {
            LinkedList<Double> cluster = clusters.get(i);
            int score = cluster.size() * 3;
            int j = i - 1;
            while (j > 0) {
                double diff = centers[i] - centers[j];
                if (diff > Consts.WATER + 0.1) {
                    break;
                }
                if (Math.abs(diff - Consts.WATER) < 0.1) {
                    score += clusters.get(j).size();
                }
                if (Math.abs(diff - 1) < 0.1) {
                    score += clusters.get(j).size();
                }

                j--;
            }
            j = i + 1;
            while (j < clusters.size()) {
                double diff = centers[j] - centers[i];
                if (diff > 1 + 0.1) {
                    break;
                }
                if (Math.abs(diff - 1) < 0.1) {
                    score += clusters.get(j).size();
                }

                j++;
            }
            if (score > bestScore) {
                bestScore = score;
                bestClusterId = i;
            }
        }


        return centers[bestClusterId];
    }

    private static List<Double> getDiffs(double[] masses, String sequence) {
        List<Double> diffsB = new ArrayList<Double>();

        double pMass = 0;
        for (int i = 0; i < sequence.length(); i++) {
            Acid acid = Acid.getAcid(sequence.charAt(i));
            if (acid == null) {
                continue;
            }
            pMass += acid.getMass();
            int cur = 0;
            while (masses[cur] - pMass < 300) {
                diffsB.add(masses[cur] - pMass);
                cur++;
                if (cur == masses.length) {
                    break;
                }
            }
        }

        Collections.sort(diffsB);
        return diffsB;
    }
}
