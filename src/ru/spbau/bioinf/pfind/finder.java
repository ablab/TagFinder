package ru.spbau.bioinf.pfind;


import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.Protein;
import ru.spbau.bioinf.tagfinder.Scan;

public class Finder {

    private static List<Protein> proteins;
    private static Map<Integer, Scan> scans;

    private static Map<Integer, Integer> ans = new HashMap<Integer, Integer>();
    private static final Set<Integer> unmatchedScans = new HashSet<Integer>();
    private static final Set<Integer> discoveredProteins = new HashSet<Integer>();

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        proteins = conf.getProteins();
        scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        
        List<MatchCandidate> candidates = new ArrayList<MatchCandidate>();

        for (Scan scan : scans.values()) {
            if (scan.getPeaks().size() > 0) {
                scan.getTags(conf);
                for (Protein protein : proteins) {
                    MatchCandidate candidate = new MatchCandidate(scan, protein);
                    if (candidate.isValid()) {
                        candidates.add(candidate);
                    }
                }
            }
        }

        /*
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        for (Map.Entry<Integer, Integer> entry : msAlignResults.entrySet()) {
            candidates.add(new MatchCandidate(scans.get(entry.getKey()), proteins.get(entry.getValue())));
        } */

        Collections.sort(candidates);
        for (int i = 0; i < 1000; i++) {
            MatchCandidate candidate = candidates.get(i);
            System.out.println(candidate.getScan().getId() + " " + candidate.getProtein().getProteinId() +  " " + candidate.getMaxScore() + " " + candidate.getScores().get(Scoring.TAG_SCORE) + " " + Math.sqrt(candidate.getScan().getPeaks().size()));
        }

        for (MatchCandidate candidate : candidates) {
            getEValueWrapper(candidate);
        }
        System.out.println("PrSM found " + ans.keySet().size());
    }


    private static int goodRequest = 0;
    private static int badRequest = 0;

    private static void getEValueWrapper(MatchCandidate candidate) throws Exception {
        int scanId = candidate.getScan().getId();
        if (ans.containsKey(scanId)) {
            return;
        }
        int proteinId = candidate.getProtein().getProteinId();
        double ret = getEValue(scanId, proteinId);
        if (ret < Configuration.EVALUE_LIMIT) {
            ans.put(scanId, proteinId);
            unmatchedScans.remove(scanId);
            discoveredProteins.add(proteinId);
            goodRequest++;
        } else {
            badRequest++;
        }
        System.out.println("<a href=http://127.0.0.1:8080/align?scanId=" + scanId + "&proteinId=" + proteinId + ">" + scanId + " " + proteinId + " " + ret + " " + goodRequest + " " + badRequest + "</a><br/>");
        if (badRequest > goodRequest) {
            throw new Exception(badRequest + " > " + goodRequest);
        }
    }

    public static double getEValue(int scanId, int proteinId) throws Exception {
            URL server = new URL("http://127.0.0.1:8080/evalue?scanId=" + scanId + "&proteinId="+proteinId);
            URLConnection conn = server.openConnection();
            BufferedReader in = new BufferedReader(
                    new InputStreamReader(
                            conn.getInputStream()));
            try {
                return Double.parseDouble(in.readLine());
            } finally {
                in.close();
            }
    }

}
