package ru.spbau.bioinf.mzscore;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import ru.spbau.bioinf.evalue.EValueServer;
import ru.spbau.bioinf.mzpeak.MzReader;
import ru.spbau.bioinf.mzpeak.MzScan;
import ru.spbau.bioinf.palign.EMassAdapter;
import ru.spbau.bioinf.tagfinder.Acid;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.Protein;
import ru.spbau.bioinf.tagfinder.Scan;

public class MzScore {

    private Protein protein;
    private List<List<double[]>> emass = new ArrayList<List<double[]>>();
    private List<Double> w = new ArrayList<Double>();

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();        
        Map<Integer,Scan> scans = conf.getScans();
        Map<Integer, MzScan> mzScans = MzReader.getMzScans();
        MzScore score = new MzScore(proteins.get(3949));
        EValueServer.init(args);
        final Map<Integer, Double> scores = new HashMap<Integer, Double>();
        for (Scan scan : scans.values()) {
            int scanId = scan.getId();
            scores.put(scanId, score.getScore(mzScans.get(scanId), scan));
        }
        
        List<Scan> res = new ArrayList<Scan>();
        res.addAll(scans.values());
        Collections.sort(res, new Comparator<Scan>() {
            public int compare(Scan o1, Scan o2) {
                double v1 = scores.get(o1.getId());
                double v2 = scores.get(o2.getId());
                double d = v1 - v2;
                if (d > 0) {
                    return -1;
                }
                if (d < 0) {
                    return 1;
                }          
                return 0;
            }
        });

        int id = 0;

        for (Scan scan : res) {
            int scanId = scan.getId();
            double v = scores.get(scanId);
            if (v == 0) {
                break;
            }

            System.out.println((++id) +  " " + scanId + " " + v + " " + EValueServer.getEvalue(scan.getId(), 3949) + " " + pos(scanId));
            if (id > 130) {
                break;
            }
        }
    }
    
    static int bad = 0;
    private static int pos(int scanId) {
        for (int i : official) {
            if (i == scanId) {
                return 0;
            }
        }
        bad++;
        return bad;
    }
    
    public MzScore(Protein protein) {
        this.protein = protein;
        String seq = protein.getSimplifiedAcids();
        String reverse = new StringBuilder(seq).reverse().toString();
        int[] atomCount = new int[5];
        atomCount[1] += 2;
        atomCount[3] += 1;

        processSequence(reverse, atomCount, 1);
        atomCount = new int[5];
        processSequence(seq, atomCount, 1);
        //atomCount = new int[5];
        //processSequence(seq.substring(1), atomCount, 1);

    }

    private void processSequence(String seq, int[] atomCount, double weight) {
        for (int i = 0; i < seq.length(); i++) {
            Acid acid = Acid.getAcid(seq.charAt(i));
            int[] add = new int[0];
            if (acid != null) {
                add = acid.getAtomCount();
            }
            for (int j = 0; j < add.length; j++) {
                atomCount[j] += add[j];
            }
            emass.add(EMassAdapter.getPeaks(atomCount));
            w.add(weight);
        }
    }


    public double getScore(MzScan mzScan, Scan scan) {
        double precursorMass = scan.getPrecursorMass();
        if (protein.getMass() + 200 < precursorMass) {
            return 0;
        }

        List<double[]> peaks = mzScan.getPeaks(precursorMass);
        for (int i = 0; i < peaks.size(); i++) {
            double[] first = peaks.get(i);
            double start = first[0];
            for (int j = i+1; j < peaks.size(); j++) {
                double[] next =  peaks.get(j);
                double v2 = (next[0] - start) * 2;
                double v3 = (next[0] - start) * 3;
                if (Math.abs(1 - v2)  < 0.04) {
                    first[1] /= 2;
                    next[1] /= 2;
                }
            }
        }

        double ans = 0;
        for (int i = 0; i < emass.size(); i++) {
            List<double[]> tpeaks = emass.get(i);
            ans += getPeptideScore(peaks, tpeaks) * w.get(i);
        }
         
        return ans;///peaks.size();
    }

    private double getPeptideScore(List<double[]> peaks, List<double[]> tpeaks) {
        double ans = 0;
        int p1 = 0;
        int p2 = 0;
        int n = tpeaks.size();
        List<Integer>[] matches = new List[n];
        for (int i = 0; i < matches.length; i++) {
            matches[i] = new ArrayList<Integer>();
        }

        while (p1 < peaks.size() && p2 < n) {
            double v1 = peaks.get(p1)[0];
            double v2 = tpeaks.get(p2)[0];
            if (Math.abs(v1-v2) < 0.1) {
                if (tpeaks.get(p2)[0] > 1) {
                    matches[p2].add(p1);
                }
                p1++;
            } else if (v1 < v2) {
                p1++;
            } else {
                p2++;
            }
        }

        double[] min = new double[n];
        double[] max = new double[n];
        for (int i = 0; i < min.length; i++) {
            double a = 10e100;
            double b = -1;
            for (int p : matches[i]) {
                double v = peaks.get(p)[1];
                if (v < a) {
                    a = v;
                }
                if (v > b) {
                    b = a;
                }
                min[i] = a;
                max[i] = b;
            }
        }
        double failed = 0;
        for (int i = 0; i < tpeaks.size(); i++) {
            double[] tp =  tpeaks.get(i);
            if (tp[1] < 25) {
                //continue;
            }
            if (matches[i].size() == 0) {
                failed += tp[1];
                continue;
            }

            if (i > 0 && i < n - 1) {
                if (matches[i-1].size() > 0 && matches[i+1].size() > 0 ) {
                    double d = 0.45 * (min[i-1] + min[i+1]) - max[i];
                    if (d > 0) {
                        failed += tp[1]/2;
                    }
                }
            }
            /*
            if (i > 0) {
                if (matches[i-1].size() > 0 && matches[i].size() > 0 ) {
                    double mid = (tpeaks.get(i - 1)[0] + tpeaks.get(i)[0])/2;
                    for (double[] peak : peaks) {                        
                        if (Math.abs(peak[0] - mid) < 0.05) {
                            failed += tp[1];
                        }
                        
                    }                    
                }
            } */
        }
        if (failed < 100) {
            ans += 100 - failed;
        }
        return ans;
    }

    private static int[] official = new int[] {
              798,
             1064,
             1066,
             1067,
             1071,
             1072,
             1076,
             1078,
             1079,
             1080,
             1082,
             1083,
             1084,
             1086,
             1087,
             1088,
             1090,
             1091,
             1092,
             1094,
             1095,
             1096,
             1098,
             1099,
             1100,
             1102,
             1103,
             1104,
             1106,
             1107,
             1108,
             1110,
             1111,
             1112,
             1114,
             1115,
             1116,
             1118,
             1119,
             1120,
             1122,
             1123,
             1126,
             1127,
             1128,
             1130,
             1131,
             1132,
             1134,
             1135,
             1136,
             1138,
             1139,
             1140,
             1142,
             1143,
             1144,
             1146,
             1147,
             1148,
             1150,
             1151,
             1152,
             1154,
             1155,
             1156,
             1158,
             1159,
             1160,
             1162,
             1163,
             1164,
             1166,
             1167,
             1168,
             1170,
             1171,
             1172,
             1174,
             1175,
             1176,
             1178,
             1179,
             1180,
             1182,
             1183,
             1184,
             1186,
             1187,
             1191,
             1194,
             1195,
             1198,
             1199,
             1200,
             1204,
             1207,
             1208,
             1214,
             1218,
             1219,
             1220,
             1222,
             1223,
             1224,
             1227,
             1228,
             1230,
             1231,
             1232,

    };    
}
