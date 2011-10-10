package ru.spbau.bioinf.tagfinder.ui;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JComponent;
import ru.spbau.bioinf.tagfinder.Acid;
import ru.spbau.bioinf.tagfinder.Analyzer;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.Peak;
import ru.spbau.bioinf.tagfinder.Scan;

public class ScanView extends JComponent {

    private static NumberFormat df = NumberFormat.getInstance();
    public static final int LINE_HEIGHT = 20;

    static {
        df.setMaximumFractionDigits(2);
    }

    List<TooltipCandidate> tooltips = new ArrayList<TooltipCandidate>();

    private Configuration conf;

    private Scan scan;
    private Dimension dimension = new Dimension(1000, 10);
    private List<List<Peak[]>> components = new ArrayList<List<Peak[]>>();

    public ScanView(Configuration conf) {
        this.conf = conf;
        addMouseMotionListener(new MouseMotionAdapter() {
            @Override
            public void mouseMoved(MouseEvent e) {
                int x = e.getX();
                int y = e.getY();
                for (TooltipCandidate tooltipCandidate : tooltips) {
                    if (tooltipCandidate.isValid(x, y)) {
                        setToolTipText(tooltipCandidate.getText());
                        break;
                    }
                }
            }
        });
    }

    public void setScan(Scan scan) {
        this.scan = scan;
        Analyzer analyzer = new Analyzer(conf);
        List<List<Peak>> components = analyzer.getComponents(scan);
        this.components.clear();
        int totalTags = 0;
        for (List<Peak> component : components) {
            List<Peak[]> tags = analyzer.getTags(component);
            totalTags += tags.size();
            this.components.add(tags);
        }

        dimension = new Dimension((int)scan.getPrecursorMass() + 200, (components.size() + totalTags) * LINE_HEIGHT + 50);
        invalidate();
    }

    @Override
    public Dimension getMinimumSize() {
        return dimension;
    }

    @Override
    public Dimension getPreferredSize() {
        return dimension;
    }

    @Override
    public Dimension getMaximumSize() {
        return dimension;
    }

    @Override
    public void paint(Graphics g) {
        tooltips.clear();
        if (scan != null) {
            List<Peak> peaks = scan.getPeaks();
            for (Peak peak : peaks) {
                int value = (int)peak.getValue();
                g.drawLine(value, 5, value, 15);
            }
            double precursorMass = scan.getPrecursorMass();
            int total = (int) precursorMass;
            g.drawLine(total, 5, total, 15);
            g.drawString(df.format(precursorMass), total + 3, 15);
        }

        int start = LINE_HEIGHT;
        for (int i = 0; i < components.size(); i++) {
            List<Peak[]> component =  components.get(i);
            double min = scan.getPrecursorMass();
            for (Peak[] tag : component) {
                double v = tag[0].getValue();
                if (v < min) {
                    min = v;
                }
            }
            g.drawString("Component " + (i + 1), 3, start + LINE_HEIGHT);
            start += LINE_HEIGHT;
            for (Peak[] tag : component) {
                for (int j = 0; j < tag.length; j++) {
                    Peak peak = tag[j];
                    int value = (int)(peak.getValue() - min);
                    g.drawLine(value, start + 3, value, start + LINE_HEIGHT);
                    tooltips.add(new TooltipCandidate(value - 5, value + 5, start, start + LINE_HEIGHT, df.format(peak.getValue())));
                    if (j + 1 < tag.length) {
                        double delta = tag[j + 1].getValue() - peak.getValue();
                        g.drawString(Acid.getAcid(delta).name(), (int)(value + delta/2 - 3), start + LINE_HEIGHT - 4);
                    }
                }
                start += LINE_HEIGHT;
            }
        }
    }
}
