package ru.spbau.bioinf.palign;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.jdom.Element;
import ru.spbau.bioinf.tagfinder.util.XmlUtil;

public class PeptideSupport {
    
    private int endPosition;

    public PeptideSupport(int endPosition) {
        this.endPosition = endPosition;
    }

    private List<PeptideModification> modifcations = new ArrayList<PeptideModification>();

    public void addModification(PeptideModification modification) {
        modifcations.add(modification);
    }

    public List<PeptideModification> getModifcations() {
        return modifcations;
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
        return xml;
    }
}
