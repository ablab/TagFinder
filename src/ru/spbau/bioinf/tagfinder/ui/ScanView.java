package ru.spbau.bioinf.tagfinder.ui;

import java.awt.Dimension;
import java.awt.Graphics;
import java.util.List;
import javax.swing.JComponent;
import ru.spbau.bioinf.tagfinder.Peak;
import ru.spbau.bioinf.tagfinder.Scan;

public class ScanView extends JComponent {

    private Scan scan;
    private Dimension dimension = new Dimension(1000, 10);

    public void setScan(Scan scan) {
        this.scan = scan;
        dimension = new Dimension((int)scan.getPrecursorMass() + 10, 100);
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
        if (scan != null) {
            List<Peak> peaks = scan.getPeaks();
            for (Peak peak : peaks) {
                int value = (int)Math.round(peak.getValue());
                g.drawLine(value, 5, value, 15);
            }
        }
    }
}
