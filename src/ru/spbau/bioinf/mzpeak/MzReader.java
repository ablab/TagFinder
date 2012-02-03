package ru.spbau.bioinf.mzpeak;


import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import javax.xml.datatype.DatatypeFactory;
import javax.xml.datatype.Duration;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import org.apache.axis.encoding.Base64;
import org.apache.log4j.Logger;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

public class MzReader {

    private static Logger log = Logger.getLogger(MzReader.class);

    private int peaksCount = 0;
    private StringBuilder charBuffer = new StringBuilder();
    private boolean compressFlag = false;
    private DefaultHandler handler = new MzXMLHandler();
    private String precision;
    private MzScan buildingScan;
    private List<MzScan> scans = new ArrayList<MzScan>();
    DatatypeFactory dataTypeFactory;


    public static Map<Integer, MzScan> getMzScans() {
        MzReader mzReader = new MzReader();
        File mzFile = new File("2006_08_16_ALS_C4_40_lipo_not_digested_30000.mzXML");
        Map<Integer, MzScan> mzScans = new HashMap<Integer, MzScan>();
        List<MzScan> mzScansList = mzReader.process(mzFile);
        for (MzScan mzScan : mzScansList) {
            mzScans.put(mzScan.getScanId(), mzScan);
        }
        return mzScans;
    }

    public static void main(String[] args) {
        File mzFile = new File("2006_08_16_ALS_C4_40_lipo_not_digested_30000.mzXML");
        new MzReader().process(mzFile);
    }

    public List<MzScan> process(File mzFile) {
        SAXParserFactory factory = SAXParserFactory.newInstance();

        try {
            dataTypeFactory = DatatypeFactory.newInstance();

            SAXParser saxParser = factory.newSAXParser();
            saxParser.parse(mzFile, handler);
        } catch (Throwable e) {
            log.error("Error reading mzXml file", e);
        }

        for (MzScan scan : scans) {
            scan.sortPoints();
        }
        return scans;
    }


    private class MzXMLHandler extends DefaultHandler {
        public void startElement(String namespaceURI, String lName, String qName, Attributes attrs) throws SAXException {
            charBuffer.setLength(0);

            if (qName.equalsIgnoreCase("scan")) {
                buildingScan = null;
                int scanNumber = Integer.parseInt(attrs.getValue("num"));
                int msLevel = Integer.parseInt(attrs.getValue("msLevel"));
                peaksCount = Integer.parseInt(attrs.getValue("peaksCount"));

                // Parse retention time
                double retentionTime = 0;
                String retentionTimeStr = attrs.getValue("retentionTime");
                if (retentionTimeStr != null) {
                    Date currentDate = new Date();
                    Duration dur = dataTypeFactory
                            .newDuration(retentionTimeStr);
                    retentionTime = dur.getTimeInMillis(currentDate) / 1000d / 60d;
                } else {
                    log.error("Failed to get retention time");
                }

                buildingScan = new MzScan(scanNumber);
                scans.add(buildingScan);
            }

            // <peaks>
            if (qName.equalsIgnoreCase("peaks")) {
                compressFlag = false;
                String compressionType = attrs.getValue("compressionType");
                if ((compressionType == null)
                        || (compressionType.equals("none")))
                    compressFlag = false;
                else
                    compressFlag = true;
                precision = attrs.getValue("precision");
            }
        }

        public void endElement(String namespaceURI, String sName, String qName) throws SAXException {
            if (qName.equalsIgnoreCase("peaks")) {

                byte[] peakBytes = Base64.decode(charBuffer.toString());

                if (compressFlag) {
                    try {
                        peakBytes = CompressionUtils.decompress(peakBytes);
                    } catch (DataFormatException e) {
                        log.error("Decompressing error", e);
                        throw new SAXException("Parsing Cancelled");
                    }
                }
                DataInputStream peakStream = new DataInputStream(
                        new ByteArrayInputStream(peakBytes));
                for (int i = 0; i < peaksCount; i++) {

                    // Always respect this order pairOrder="m/z-int"
                    double massOverCharge;
                    double intensity;
                    try {
                        if ("64".equals(precision)) {
                            massOverCharge = peakStream.readDouble();
                            intensity = peakStream.readDouble();
                        } else {
                            massOverCharge = (double) peakStream.readFloat();
                            intensity = (double) peakStream.readFloat();
                        }
                    } catch (IOException e) {
                        log.error("Error reading peaks data", e);
                        throw new SAXException("Parsing Cancelled");
                    }
                    buildingScan.addPoint(massOverCharge, intensity);
                }
        }
    }

    public void characters(char buf[], int offset, int len)
            throws SAXException {
        charBuffer.append(buf, offset, len);
    }
}


public static class CompressionUtils {

    public static byte[] decompress(byte compressedBytes[])
            throws DataFormatException {

        Inflater decompresser = new Inflater();

        decompresser.setInput(compressedBytes);

        byte[] resultBuffer = new byte[compressedBytes.length * 2];
        byte[] resultTotal = new byte[0];

        int resultLength = decompresser.inflate(resultBuffer);

        while (resultLength > 0) {
            byte previousResult[] = resultTotal;
            resultTotal = new byte[resultTotal.length + resultLength];
            System.arraycopy(previousResult, 0, resultTotal, 0,
                    previousResult.length);
            System.arraycopy(resultBuffer, 0, resultTotal,
                    previousResult.length, resultLength);
            resultLength = decompresser.inflate(resultBuffer);
        }

        decompresser.end();

        return resultTotal;
    }

}


}
