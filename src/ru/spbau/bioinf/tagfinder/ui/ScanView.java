package ru.spbau.bioinf.tagfinder.ui;

import java.awt.Color;
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
import ru.spbau.bioinf.tagfinder.PrecursorMassShiftFinder;
import ru.spbau.bioinf.tagfinder.Protein;
import ru.spbau.bioinf.tagfinder.Scan;
import ru.spbau.bioinf.tagfinder.ShiftEngine;

public class ScanView extends JComponent {

    private static NumberFormat df = NumberFormat.getInstance();
    public static final int LINE_HEIGHT = 20;

    static {
        df.setMaximumFractionDigits(2);
    }

    List<TooltipCandidate> tooltips = new ArrayList<TooltipCandidate>();

    private Configuration conf;

    private Scan scan;
    int proteinId = -1;
    private Protein protein = null;

    private Dimension dimension = new Dimension(1000, 10);
    private List<List<Peak[]>> components = new ArrayList<List<Peak[]>>();

    private List<Protein> proteins;
    private double[] proteinSpectrum;

    public ScanView(Configuration conf, List<Protein> proteins) {
        this.conf = conf;
        this.proteins = proteins;
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

    public boolean setProteinId(int newProteinId) {
        if (newProteinId != proteinId) {
            if (newProteinId >= 0 && proteins.size() > newProteinId) {
                this.proteinId = newProteinId;
                protein = proteins.get(proteinId);
                proteinSpectrum = ShiftEngine.getSpectrum(protein.getSimplifiedAcids());
                return true;
            }
        }
        return false;
    }

    public Protein getProtein() {
        return protein;
    }

    public boolean setScan(Scan scan) {
        if (scan != this.scan) {
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

            dimension = new Dimension((int)scan.getPrecursorMass() + 200, (components.size() + totalTags + 6) * LINE_HEIGHT);
            return true;
        }
        return false;
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
        List<Peak> peaks = null;
        if (scan != null) {
            peaks = scan.getPeaks();
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

        double bestShift = 0;

        if (protein != null) {
            String sequence = protein.getSimplifiedAcids();
            double pos = 0;
            for (int cur = 0; cur < sequence.length(); cur++) {
                pos += drawLetter(g, start, pos, Acid.getAcid(sequence.charAt(cur)));
                drawLine(g, start, pos);
            }
            start+= LINE_HEIGHT;

            if (scan != null) {
                double precursorMass =  scan.getPrecursorMass() + PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan);
                List<Double> shifts = ShiftEngine.getShifts(peaks, precursorMass, proteinSpectrum);

                double bestScore = 0;
                double[] spectrum = ShiftEngine.getSpectrum(peaks, precursorMass);
                for (Double shift : shifts) {
                    double nextScore = ShiftEngine.getScore(spectrum, proteinSpectrum, shift);
                    if (nextScore > bestScore) {
                        bestScore = nextScore;
                        bestShift = shift;
                    }
                }
            }
        }

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
                    double peakValue = peak.getValue();
                    double value = (peakValue - min);
                    if (proteinSpectrum != null) {
                        g.setColor(ShiftEngine.contains(proteinSpectrum, peakValue + bestShift) ? Color.GREEN : Color.BLACK);
                    }
                    drawLine(g, start, value);
                    g.setColor(Color.BLACK);
                    tooltips.add(new TooltipCandidate(value - 5, value + 5, start, start + LINE_HEIGHT, df.format(peakValue)));
                    if (j + 1 < tag.length) {
                        double delta = tag[j + 1].getValue() - peakValue;
                        drawLetter(g, start, value, Acid.getAcid(delta));
                    }
                }
                start += LINE_HEIGHT;
            }
        }
    }

    private void drawLine(Graphics g, int start, double value) {
        g.drawLine((int)value, start + 3, (int)value, start + LINE_HEIGHT);
    }

    private double drawLetter(Graphics g, int start, double pos, Acid acid) {
        double delta = acid.getMass();
        g.drawString(acid.name(), (int) (pos + delta / 2 - 3), start + LINE_HEIGHT - 4);
        return delta;
    }
}
