package ru.spbau.bioinf.palign;

import java.util.ArrayList;
import java.util.List;
import org.jdom.Element;
import ru.spbau.bioinf.tagfinder.Acid;
import ru.spbau.bioinf.tagfinder.util.XmlUtil;

public class Cleavage {
    
    private int position;
    private int direction;

    private double modification = 0;
    
    private String seq;

    public Cleavage(int position, int direction, double modification, String seq) {
        this.position = position;
        this.direction = direction;
        this.modification = modification;
        this.seq = seq;
    }

    private List<PeptideSupport> supports = new ArrayList<PeptideSupport>();

    public void addSupport(PeptideSupport support) {
        supports.add(support);
    }

    public List<PeptideSupport> getSupports() {
        return supports;
    }

    public int getPosition() {
        return position;
    }

    public int getDirection() {
        return direction;
    }

    public double getModification() {
        return modification;
    }

    public Element toXml() {
        Element xml = new Element("cleavage");
        int[] atomCount = new int[5];
        double mass = modification;
        for (int i = 0; i < seq.length(); i++) {
            Acid acid = Acid.getAcid(seq.charAt(i));
            int[] add = new int[0];
            if (acid != null) {
                mass += acid.getMass();
                add = acid.getAtomCount();
            }
            for (int j = 0; j < add.length; j++) {
                atomCount[j] += add[j];
            }
            Element peptide = new Element("peptide");
            xml.addContent(peptide);
            XmlUtil.addElement(peptide, "mass", mass);
            XmlUtil.addElement(peptide, "acid", acid.name());
            PeptideSupport.addPeaks(peptide, atomCount);
        }
        XmlUtil.addElement(xml, "position", position);
        XmlUtil.addElement(xml, "direction", direction);
        XmlUtil.addElement(xml, "modification", modification);        
        for (PeptideSupport support : supports) {
            xml.addContent(support.toXml());
        }
        return xml;
    }
}
