package ru.spbau.bioinf.palign;

import org.jdom.Element;
import ru.spbau.bioinf.mzpeak.MzPoint;
import ru.spbau.bioinf.tagfinder.util.XmlUtil;

public class MzSupport {
    
    private MzPoint mzPoint;
    private int charge;
    private double error;

    public MzSupport(MzPoint mzPoint, int charge, double error) {
        this.mzPoint = mzPoint;
        this.charge = charge;
        this.error = error;
    }

    public Element toXml() {
        Element xml = new Element("mz-support");
        XmlUtil.addElement(xml, "mass", mzPoint.getMass());
        XmlUtil.addElement(xml, "intencity", mzPoint.getIntencity());
        XmlUtil.addElement(xml, "charge", charge);
        XmlUtil.addElement(xml, "error", error);
        return xml;
    }

}
