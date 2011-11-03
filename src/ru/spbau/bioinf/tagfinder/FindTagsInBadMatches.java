package ru.spbau.bioinf.tagfinder;

import java.io.File;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class FindTagsInBadMatches {

    private static NumberFormat df = NumberFormat.getInstance();
    private static List<Protein> proteins;
    private static Configuration conf;

    static {
        df.setMaximumFractionDigits(2);
    }

    public static void main(String[] args) throws Exception {
        conf = new Configuration(args);
        proteins = conf.getProteins();
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        //Map<Integer, Integer> badMsAlignResults = conf.getBadMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        for (int i = 5; i < 12; i++) {
            PrintWriter out = ReaderUtil.createOutputFile(new File("tags", "tags_" + i + "_more.txt"));
            process(out, i, msAlignResults, scans);
            out.close();
        }
    }

    private static void process(PrintWriter out, int len, Map<Integer, Integer> msAlignResults, Map<Integer, Scan> scans) {
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);

        for (Scan scan : scans.values()) {
            int scanId = scan.getId();
            if (msAlignResults.containsKey(scanId)) {
                continue;
            }
            //if (!badMsAlignResults.containsKey(scanId)) {
                //continue;
            //}
            double shift =
                    //0
                    PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan)
                    ;

            List<Peak> peaks = //msAlignPeaks.get(scanId);
                                    scan.createSpectrumWithYPeaks(shift);
                                    //scan.createStandardSpectrum();

                //List<Peak> yPeaks = new ArrayList<Peak>();
                Collections.sort(peaks);
//                for (Peak peak : peaks) {
//                    yPeaks.add(peak.getYPeak(scan.getPrecursorMass()));
//                }
//                Collections.sort(yPeaks);
            process(out, len, scan, peaks);
            //System.out.println("Reverse tags: ");
            //process(scan, yPeaks, 1);
        }
    }

    private static void process(PrintWriter out, int len, Scan scan, List<Peak> peaks) {
        //System.out.println("scan = " + scan.getId());
        GraphUtil.generateGapEdges(conf, peaks, 1);
        used.clear();
        for (Peak peak : peaks) {
            checkTags(len, scan, peak, peak.getMass(), "");
        }
        for (Protein protein : proteins) {
            int score = 0;
            String text ="";
            for (String s : used) {
                if (protein.getSimplifiedAcids().contains(s)) {
                    score ++;
                    text += s + " ";
                }
            }
            if (score > 0) {
                out.println("Match " + score + " " + text +  scan.getId() + " " + protein.getProteinId());
            }
        }
    }
    static Set<String> used = new HashSet<String>();
    
     public static void checkTags(int len, Scan scan, Peak peak, double startMass, String prefix) {
        if(prefix.length() == len) {
            if (used.contains(prefix)) {
                return;
            } else {
                used.add(prefix);
            }

            String text = scan.getId() + " " + prefix;
            boolean add = false;

            for (Protein protein : proteins) {
                int cur = -1;
                String sequence = protein.getSimplifiedAcids();
                while ((cur = sequence.indexOf(prefix, cur + 1))>=0) {
                    double mass = 0;
                    for (int i = 0; i < cur - 1; i++) {
                        try {
                            mass += Acid.getAcid(sequence.charAt(i)).getMass();
                        } catch (Throwable e) {
                            //System.out.println("Wrong symbol in protein " + protein.getProteinId());
                        }
                    }
                    add = true;
                    text += " " + protein.getProteinId();
                    
                    //if (Math.abs(mass - startMass) < scan.getPrecursorMass() * TEN_PPM * 1.5) {
                        //System.out.println(scan.getId() + " " + protein.getProteinId() + " " + prefix + " " + (mass - startMass));
                    //}
                }
            }
            if (add) {
                //System.out.println(text);
            }
            return;
        }

        for (Peak next : peak.getNext()) {
            for (Acid acid : Acid.values()) {
                if (acid.match(conf.getEdgeLimits(peak, next))) {
                    checkTags(len, scan, next, startMass, prefix + acid.name());
                }
            }
        }
    }
}
