package ru.spbau.bioinf.tagfinder.view;

import org.jdom.Element;
import ru.spbau.bioinf.tagfinder.Acid;
import ru.spbau.bioinf.tagfinder.Peak;

import java.util.HashMap;
import java.util.Set;

public class Table {

    private HashMap<Integer, Row> rows = new HashMap<Integer, Row>();

    public static int[] getRange(HashMap<Integer, ?> map) {
        int[] ans = createInitialRange();
        return getRange(map, ans);
    }

    private static int[] createInitialRange() {
        return new int[]{Integer.MAX_VALUE, Integer.MIN_VALUE};
    }

    public static int[] getRange(HashMap<Integer, ?> map, int[] ans) {
        Set<Integer> values = map.keySet();
        return getRange(values, ans);
    }

    private static int[] getRange(Set<Integer> values, int[] ans) {
        for (int v : values) {
            if (ans[0] > v) {
                ans[0] = v;
            }
            if (ans[1] < v) {
                ans[1] = v;
            }
        }
        return ans;
    }

    public void addTag(int row, int col, Peak[] peaks) {
        if (!rows.containsKey(row)) {
            rows.put(row, new Row(row));
        }
        Row r = rows.get(row);
        for (int i = 0; i < peaks.length - 1; i++) {
            Peak peak = peaks[i];
            r.addCell(col, new PeakContent(peak));
            col++;
            Acid acid = Acid.getAcid(peaks[i+1].diff(peak));
            r.addCell(col, new AcidContent(acid));
            col++;
        }
        r.addCell(col, new PeakContent(peaks[peaks.length - 1]));
    }

    public Element toXML() {
        Element table = new Element("table");
        int[] rLimits = getRange(rows);
        int[] cLimits = createInitialRange();
        for (Row row : rows.values()) {
            cLimits = getRange(row.getCells(), cLimits);
        }
        for (int r = rLimits[0]; r <= rLimits[1]; r++) {
            Element row = new Element("row");
            for (int c = cLimits[0]; c <= cLimits[1];  c++) {
                Element cell = new Element("cell");
                if (rows.containsKey(r)) {
                    Content content = rows.get(r).getContent(c);
                    if (content != null) {
                        cell.addContent(content.toXml());
                    }
                }
                row.addContent(cell);
            }
            table.addContent(row);
        }
        return table;
    }
}
