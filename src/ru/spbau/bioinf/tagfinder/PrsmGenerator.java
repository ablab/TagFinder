package ru.spbau.bioinf.tagfinder;

import java.io.File;
import java.util.List;
import java.util.Map;
import org.jdom.Document;
import org.jdom.Element;
import ru.spbau.bioinf.tagfinder.util.XmlUtil;

public class PrsmGenerator {
    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        Map<Integer,Scan> scans = conf.getScans();
        List<Protein> proteins = conf.getProteins();
        Map<Integer, Integer> results = conf.getMSAlignResults();
        for (int scanId : scans.keySet()) {
            Document doc = new Document();
            Element prsm = new Element("prsm");
            doc.setRootElement(prsm);
            prsm.addContent(scans.get(scanId).toXml());
            Integer proteinId = results.get(scanId);
            if (proteinId != null) {
                prsm.addContent(proteins.get(proteinId).toXml());
            }
            XmlUtil.saveXml(doc, new File(conf.getXmlPrsmDir(),"scan" + scanId + ".xml"));
        }
    }
}
