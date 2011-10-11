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
import ru.spbau.bioinf.tagfinder.Consts;
import ru.spbau.bioinf.tagfinder.Peak;
import ru.spbau.bioinf.tagfinder.PeakType;
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
    double bestShift = 0;

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
                if (scan != null) {
                    initBestShift();
                } else {
                    bestShift = 0;
                }
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

            if (proteinSpectrum != null) {
                initBestShift();
            }

            dimension = new Dimension((int)scan.getPrecursorMass() + 400, (components.size() + totalTags + 6) * LINE_HEIGHT);
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
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, getWidth(), getHeight());
        g.setColor(Color.BLACK);
        tooltips.clear();
        int y = 0;
        List<Peak> peaks = null;
        if (scan != null) {
            peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan));
            for (Peak peak : peaks) {
                drawPeak(g, peak, y, -bestShift);
            }
            double precursorMass = scan.getPrecursorMass();
            double end = precursorMass + bestShift;
            drawPeak(g, y, bestShift);
            drawPeak(g, y, end);
            g.drawString(df.format(precursorMass) +  " + " + df.format(bestShift), (int)end + 3, 15);
        }

        y += LINE_HEIGHT;



        if (protein != null) {
            String sequence = protein.getSimplifiedAcids();
            double pos = 0;
            for (int cur = 0; cur < sequence.length(); cur++) {
                pos += drawLetter(g, y, pos, Acid.getAcid(sequence.charAt(cur)));
                drawPeak(g, y, pos);
            }
            y+= LINE_HEIGHT;
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
            g.drawString("Component " + (i + 1), 3, y + LINE_HEIGHT);
            y += LINE_HEIGHT;
            for (Peak[] tag : component) {
                for (int j = 0; j < tag.length; j++) {
                    Peak peak = tag[j];
                    double peakValue = drawPeak(g, peak, y, min);
                    if (j + 1 < tag.length) {
                        double delta = tag[j + 1].getValue() - peakValue;
                        drawLetter(g, y, peakValue - min, Acid.getAcid(delta));
                    }
                }
                y += LINE_HEIGHT;
            }
        }
    }

    private double drawPeak(Graphics g, Peak peak, int y, double min) {
        double peakValue = peak.getValue();
        double value = (peakValue - min);
        if (proteinSpectrum != null) {
            Color color = Color.BLACK;
            double v = peakValue + bestShift;
            if (ShiftEngine.contains(proteinSpectrum, v)) {
                color = Color.GREEN;
            } else if (ShiftEngine.contains(proteinSpectrum, v, v - 1, v + 1, v + Consts.WATER, v - Consts.WATER, v - Consts.AMMONIA, v + Consts.AMMONIA)) {
                color = Color.ORANGE;
            }
            g.setColor(color);
        }
        drawPeak(g, y, value, peakValue);
        drawIon(g, y, value, peak.getPeakType());
        g.setColor(Color.BLACK);
        return peakValue;
    }

    private void initBestShift() {
        double precursorMass =  scan.getPrecursorMass() + PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan);
        List<Double> shifts = ShiftEngine.getShifts(scan.getPeaks(), precursorMass, proteinSpectrum);

        double bestScore = 0;
        double[] spectrum = ShiftEngine.getSpectrum(scan.getPeaks(), precursorMass);
        for (Double shift : shifts) {
            double nextScore = ShiftEngine.getScore(spectrum, proteinSpectrum, shift);
            if (nextScore > bestScore) {
                bestScore = nextScore;
                bestShift = shift;
            }
        }
    }

    private void drawPeak(Graphics g, int y, double value) {
        drawPeak(g, y, value, value);
    }

    private void drawPeak(Graphics g, int y, double value, double peakValue) {
        g.drawLine((int)value, y + 3, (int)value, y + LINE_HEIGHT);
        tooltips.add(new TooltipCandidate(value - 5, value + 5, y, y + LINE_HEIGHT, df.format(peakValue)));
    }

    private void drawIon(Graphics g, int y, double value, PeakType peakType) {
        int v = (int) value;
        int x1 = peakType == PeakType.B ? v - 3 : v + 3;
        g.drawLine(x1, y + LINE_HEIGHT, v, y + LINE_HEIGHT);
        g.drawLine(x1, y + 3, v, y + 3);
    }

    private double drawLetter(Graphics g, int start, double pos, Acid acid) {
        double delta = acid.getMass();
        g.drawString(acid.name(), (int) (pos + delta / 2 - 3), start + LINE_HEIGHT - 4);
        return delta;
    }
}
