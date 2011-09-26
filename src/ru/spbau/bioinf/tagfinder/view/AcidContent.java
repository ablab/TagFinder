package ru.spbau.bioinf.tagfinder.view;

import org.jdom.Element;
import ru.spbau.bioinf.tagfinder.Acid;
import ru.spbau.bioinf.tagfinder.util.XmlUtil;

public class AcidContent implements Content{

    private Acid acid;

    public AcidContent(Acid acid) {
        this.acid = acid;
    }


    @Override
    public Element toXml() {
        Element ans = new Element("acid");
        XmlUtil.addElement(ans, "name", acid.name());
        return ans;

    }
}
