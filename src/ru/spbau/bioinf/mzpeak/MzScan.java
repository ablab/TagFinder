package ru.spbau.bioinf.mzpeak;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class MzScan {
    
    private int scanId;

    private List<MzPoint> points = new ArrayList<MzPoint>();

    public MzScan(int scanNamber) {
        this.scanId = scanNamber;
    }

    public int getScanId() {
        return scanId;
    }

    public void sortPoints() {
        Collections.sort(points);
    }

    public void addPoint(double mass, double intencity) {
        points.add(new MzPoint(mass, intencity));
    }

    public List<MzPoint> getPoints() {
        return points;
    }
    
    public List<Double> getSupportedPoints(double ppmError, int dir) {
        double[] a = new double[30 * points.size()];
        for (int i = 0; i < points.size(); i++) {
            MzPoint mzPoint = points.get(i);
            double mass = mzPoint.getMass();
            for (int j = 1; j <= 30; j++) {
                a[i*30 + j - 1] = mass * j - j ;
            }
        }
        Arrays.sort(a);
        int ans = 0;
        List<Double> masses = new ArrayList<Double>();
        List<Double> centers = new ArrayList<Double>();
        for (int i = 0; i < a.length - 1; i++) {
            if (a[i + 1] -  a[i] < ppmError * (a[i+1] + a[i])) {
                if (masses.isEmpty()) {
                    masses.add(a[i]);
                }
                masses.add(a[i + 1]);
            } else {
                if (masses.size() > 0) {
                    centers.add(getCenter(masses));
                    masses.clear();
                }
            }
            
        }

        for (Iterator<Double> iterator = centers.iterator(); iterator.hasNext(); ) {
            Double center = iterator.next();
            boolean supported = false;
            for (int i = 0; i < a.length; i++) {
                if (Math.abs(Math.abs(a[i] - center) - 1) < ppmError * 2 * center) {
                    supported = true;
                    break;
                }
            }
            if (!supported) {
                iterator.remove();
            }
        }
        
        List<Double> extend = new ArrayList<Double>();
        for (double center : centers) {
            double pos = center - 20;
            int cur = 0;
            while (pos < center + 20 && cur < a.length) {
                double diff = pos - a[cur];
                if (Math.abs(diff) < ppmError * (pos + a[cur])) {
                    extend.add(pos);
                    pos += 1;
                    cur++;
                } else {
                    if (diff > 0) {
                        cur++;
                    } else {
                        pos +=1;
                    }
                }
            }
            
        }
        
        //extend = filterScale(extend, ppmError);
        return extend;
    }

    public static List<Double> filterScale(List<Double> masses, double ppmError) {
        Collections.sort(masses);
        List<Double> ans = new ArrayList<Double>();

        for (int i = 0; i < masses.size(); i++) {
            double m1 = masses.get(i);
            boolean correct = true;
            for (int j = i + 1; j < masses.size(); j++) {
                double m2 =  masses.get(j);
                for (int k = 2; k <= 30 ; k++) {
                     if (Math.abs(k * m1 - m2) < ppmError * 2 * m2) {
                         correct = false;
                         break;
                     }
                }
                if (!correct) {
                    break;
                }
            }
            if (correct) {
                ans.add(m1);
            }
        }
        return ans;
    }
    
    public static double getCenter(List<Double> masses) {
        double total = 0;
        for (Double mass : masses) {
            total += mass;
        }
        return total / masses.size();
        
    }
}
