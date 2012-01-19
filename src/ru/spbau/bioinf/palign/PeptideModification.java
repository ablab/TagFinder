package ru.spbau.bioinf.palign;

import org.jdom.Element;
import ru.spbau.bioinf.tagfinder.Peak;
import ru.spbau.bioinf.tagfinder.util.XmlUtil;

public class PeptideModification {
    
    private Peak peak;

    private ModificationType type;

    private double error;

    public PeptideModification(Peak peak, ModificationType type, double error) {
        this.peak = peak;
        this.type = type;
        this.error = error;
    }

    public ModificationType getType() {
        return type;
    }

    public Peak getPeak() {
        return peak;
    }

    public Element toXml() {
        Element xml = new Element("modification");
        xml.addContent(peak.toXml());
        XmlUtil.addElement(xml, "modification-type", type.name());
        XmlUtil.addElement(xml, "error", error);
        return xml;
    }
}
