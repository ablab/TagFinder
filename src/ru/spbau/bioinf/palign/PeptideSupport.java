package ru.spbau.bioinf.palign;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jdom.Element;
import ru.spbau.bioinf.mzpeak.MzPoint;
import ru.spbau.bioinf.mzpeak.MzScan;
import ru.spbau.bioinf.tagfinder.util.XmlUtil;

public class PeptideSupport {

    private double expectedMass;
    private int endPosition;
    private int[] atomsCount;
    private double proteinModification;
    
    private Map<Integer, List<MzSupport>> mzSupports = new HashMap<Integer, List<MzSupport>>();

    public PeptideSupport(double expectedMass, int endPosition, int[] atomsCount, double proteinModification) {
        this.expectedMass = expectedMass;
        this.endPosition = endPosition;
        this.atomsCount = atomsCount;
        this.proteinModification = proteinModification;
    }

    private List<PeptideModification> modifcations = new ArrayList<PeptideModification>();

    public void addModification(PeptideModification modification) {
        modifcations.add(modification);
    }

    public List<PeptideModification> getModifcations() {
        return modifcations;
    }

    public void addScanInfo(MzScan scan) {
        mzSupports.clear();
        List<MzPoint> points = scan.getPoints();
        mzSupports.put(0, new ArrayList<MzSupport>());
        for (MzPoint point : points) {
            double mz = point.getMass();
            for (int charge = 1; charge <= 30; charge++) {
                double mass = charge * mz - charge;
                if (Math.abs(mass - expectedMass) < 15) {
                    double error = mass - (expectedMass);
                    mzSupports.get(0).add(new MzSupport(point, charge, error));
                }
                /*
                for (int delta = -15; delta < 15; delta++) {
                    double error = mass - (expectedMass + delta);
                    if (Math.abs(error) < expectedMass * 0.0001) {
                        if (!mzSupports.containsKey(delta)) {
                            mzSupports.put(delta, new ArrayList<MzSupport>());
                        }
                        mzSupports.get(delta).add(new MzSupport(point, charge, error));
                    }
                }
                */
            }
        }
    }

    public Element toXml() {
        Element xml = new Element("support");
        XmlUtil.addElement(xml, "end-position", endPosition);
        Collections.sort(modifcations, new Comparator<PeptideModification>() {
            public int compare(PeptideModification o1, PeptideModification o2) {
                return o1.getType().compareTo(o2.getType());
            }
        });
        for (PeptideModification modifcation : modifcations) {
            xml.addContent(modifcation.toXml());
        }

        XmlUtil.addElement(xml, "expected-mass", expectedMass);
        XmlUtil.addElement(xml, "protein-modification", proteinModification);
        List<Integer> deltas = new ArrayList<Integer>();
        deltas.addAll(mzSupports.keySet());
        Collections.sort(deltas);
        for (Integer delta : deltas) {
            Element peak = new Element("peak");
            xml.addContent(peak);
            XmlUtil.addElement(peak, "delta", delta);
            for (MzSupport mzSupport : mzSupports.get(delta)) {
                peak.addContent(mzSupport.toXml());
            }
        }

        addPeaks(xml, atomsCount);

        return xml;        
    }

    public static void addPeaks(Element xml, int[] atomsCount) {
        ArrayList<double[]> peaks = EMassAdapter.getPeaks(atomsCount);
        for (double[] peak : peaks) {
            Element tpeak = new Element("tpeak");
            XmlUtil.addElement(tpeak, "value", peak[0]);
            XmlUtil.addElement(tpeak, "score", peak[1]);
            xml.addContent(tpeak);
        }
    }
}
