package ru.spbau.bioinf.tagfinder.stat;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import ru.spbau.bioinf.tagfinder.Acid;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.GraphUtil;
import ru.spbau.bioinf.tagfinder.Peak;
import ru.spbau.bioinf.tagfinder.Protein;
import ru.spbau.bioinf.tagfinder.Scan;

public class PlaceStatistics {

    public static final double EPSILON = 0.1;
    private Configuration conf;

    public PlaceStatistics(Configuration conf) {
        this.conf = conf;
    }

    private static NumberFormat df = NumberFormat.getInstance();

    static {
        df.setMaximumFractionDigits(2);
    }

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        PlaceStatistics placeStatistics = new PlaceStatistics(conf);
        Set<Integer> usedProteins = new HashSet<Integer>();
        Map<Integer, List<Peak>> msAlignPeaks = conf.getMSAlignPeaks();
        Map<Integer, Map<Double, String>> msAlignDatas = conf.getMSAlignData();
        for (int key : keys) {
            Scan scan = scans.get(key);
            int scanId = scan.getId();
            if (msAlignResults.containsKey(scanId)) {
                Integer proteinId = msAlignResults.get(scanId);
                if (usedProteins.contains(proteinId)) {
                    continue;
                }
                usedProteins.add(proteinId);
                if (scanId == 1946) {
                    continue;
                }


                String sequence = proteins.get(proteinId).getSimplifiedAcids();
                List<Peak> peaks = msAlignPeaks.get(scanId);
                GraphUtil.generateEdges(conf, peaks);
                Map<Double, String> msAlignData = msAlignDatas.get(scanId);
                int[][] stat = new int[100][3];
                Set<String> tags = GraphUtil.generateTags(conf, peaks);
                for (String tag : tags) {
                    int len = tag.length();
                    double prev = -1;
                    for (Peak peak : peaks) {
                        if (peak.getValue() - prev < 0.02) {
                            continue;
                        }
                        double pos = peak.getValue();
                        int i;
                        for (i = 0; i < tag.length(); i++) {
                            pos += Acid.getAcid(tag.charAt(i)).getMass();
                            boolean found = false;
                            for (Peak p2 : peaks) {
                                if (Math.abs(p2.getValue() - pos) < EPSILON) {
                                    found = true;
                                    break;
                                }
                            }
                            if (!found) {
                                break;
                            }
                        }
                        if (i == tag.length()) {
                            stat[len][0]++;
                            if (sequence.contains(tag)) {
                                stat[len][1]++;
                                for (Map.Entry<Double, String> entry : msAlignData.entrySet()) {
                                    if (entry.getValue().startsWith(tag) && Math.abs(peak.getValue() - entry.getKey()) < EPSILON) {
                                        stat[len][2]++;
                                        prev = peak.getValue();
                                    }
                                }
                            }
                        }
                    }
                }

                System.out.print(scanId + " " +  proteinId);
                for (int i = 1; i < stat.length; i++) {
                    long total = stat[i][0];
                    if (total > 0) {
                        System.out.print(" " + df.format((100d * stat[i][1])/total));
                        System.out.print(" " + df.format((100d * stat[i][2])/total));
                    } else {
                        break;
                    }
                }
                System.out.println();

            }
        }
    }


    public void generateTags(Set<String> tags, List<Peak> peaks) {
        for (Peak peak : peaks) {
            generateTags(tags, "", peak);
        }
    }

    public void generateTags(Set<String> tags, String prefix, Peak peak) {
        tags.add(prefix);
        for (Peak next : peak.getNext()) {
            for (Acid acid : Acid.values()) {
                if (acid.match(conf.getEdgeLimits(peak, next))) {
                    generateTags(tags, prefix + acid.name(), next);
                }
            }
        }
    }

}
