package ru.spbau.bioinf.tagfinder;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class FindUnidentifiedWithGoodTags {
    public static void main(String[] args) throws Exception {
        File cacheFile = new File("cache.txt");
        BufferedReader in = ReaderUtil.createInputReader(cacheFile);
        Set<Integer> identified = new HashSet<Integer>();
        Map<Integer, Integer> ans = new HashMap<Integer, Integer>();
        do {
            String s = in.readLine();
            if (s == null) {
                break;
            }
            String[] data = s.split(" ");
            double eValue = Double.parseDouble(data[2]);
            if (eValue < Configuration.EVALUE_LIMIT) {
                int scanId = Integer.parseInt(data[0]);
                identified.add(scanId);    
                ans.put(scanId, Integer.parseInt(data[1]));
            }            
        } while (true);
        Configuration conf = new Configuration(args);
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        List<Protein> proteins = conf.getProteins();
        int good = 0;
        int bad = 0;
        //keys.clear();
        //keys.add(1718);
        for (int scanId : keys) {
            if (!identified.contains(scanId)) {
                Scan scan = scans.get(scanId);
                List<Peak> peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan));
                Collections.sort(peaks);
                //GraphUtil.generateEdges(conf, peaks);
                GraphUtil.generateGapEdges(conf, peaks, 1);
                //ValidTags.filterMonotags(peaks);
                Peak[] tagPeaks = GraphUtil.findBestTag(peaks);
                if (tagPeaks.length > 8) {
                    String tag = "";
                    for (int j = 0; j < tagPeaks.length  - 1; j++) {
                        tag += Acid.getAcid(tagPeaks[j + 1].getValue() - tagPeaks[j].getValue()).name();
                    }
                    System.out.println(scanId + " " + (tagPeaks.length - 1) + //" " + scan.getPeaks().size() +
                            " " + tag + " " + new StringBuilder(tag).reverse().toString());
                }
            } else {
                /*
                int proteinId = ans.get(scanId);
                System.out.println(scanId + " " + proteinId);
                Protein protein = proteins.get(proteinId);
                String sequence = protein.getSimplifiedAcids();

                Scan scan = scans.get(scanId);
                List<Peak> peaks = scan.createSpectrumWithYPeaks(0);
                Collections.sort(peaks);                
                GraphUtil.generateGapEdges(conf, peaks, 1);
                Map<String, Peak> tagsMap = GraphUtil.generateTagsWithStarts(conf, peaks);
                boolean isGood = false;

                for (String tag : tagsMap.keySet()) {
                    if (tag.length() == 3) {
                        if (sequence.contains(tag)) {
                            isGood = true;
                            break;
                        }
                    }
                }
                if (isGood) {
                    good++;
                } else {
                    bad++;
                }
                  */
            }
        }
        System.out.println("bad = " + bad);
        System.out.println("good = " + good);
        System.out.println("total = " + ans.keySet().size());
        
    }    
}
