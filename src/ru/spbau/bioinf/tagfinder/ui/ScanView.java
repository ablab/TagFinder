package ru.spbau.bioinf.tagfinder.ui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.swing.JComponent;
import ru.spbau.bioinf.tagfinder.Acid;
import ru.spbau.bioinf.tagfinder.Analyzer;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.Consts;
import ru.spbau.bioinf.tagfinder.MassAlign;
import ru.spbau.bioinf.tagfinder.MassMatch;
import ru.spbau.bioinf.tagfinder.ModificationDiscoverer;
import ru.spbau.bioinf.tagfinder.Peak;
import ru.spbau.bioinf.tagfinder.PeakType;
import ru.spbau.bioinf.tagfinder.PrecursorMassShiftFinder;
import ru.spbau.bioinf.tagfinder.Protein;
import ru.spbau.bioinf.tagfinder.Scan;
import ru.spbau.bioinf.tagfinder.ShiftEngine;
import ru.spbau.bioinf.tagfinder.stat.PlaceStatistics;

public class ScanView extends JComponent {

    private static NumberFormat df = NumberFormat.getInstance();
    public static final int LINE_HEIGHT = 20;

    static {
        df.setMaximumFractionDigits(2);
        df.setGroupingUsed(false);
    }

    List<TooltipCandidate> tooltips = new ArrayList<TooltipCandidate>();

    private final Configuration conf;
    private final TagFinder tagFinder;

    private Scan scan;
    List<Peak> peaks = null;

    int proteinId = -1;
    private Protein protein = null;
    private double scale = 0.5;
    private Peak selectedPeak = null;

    private Dimension dimension = new Dimension(1000, 10);
    private List<List<Peak[]>> components = new ArrayList<List<Peak[]>>();

    private double[] proteinSpectrum;
    private double bestShiftB = 0;
    private double bestShiftY = 0;
    private List<Double> rangeShifts = new ArrayList<Double>();
    private List<MassMatch> goodMatches;
    private List<MassMatch> badMatches;
    private double bestDiff;

    public ScanView(Configuration conf, final TagFinder tagFinder) {
        this.conf = conf;
        this.tagFinder = tagFinder;

        addMouseMotionListener(new MouseAdapter() {
            @Override
            public void mouseMoved(MouseEvent e) {
                TooltipCandidate tooltip = getTooltip(e);
                if (tooltip != null) {
                    setToolTipText(tooltip.getText());
                }
            }
        });

        addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                requestFocus();
                TooltipCandidate tooltip = getTooltip(e);
                if (tooltip != null) {
                    if (tooltip.getPeak() != null) {
                        if (selectedPeak != tooltip.getPeak()) {
                            selectedPeak = tooltip.getPeak();
                            tagFinder.updateStatus("Mass: " + selectedPeak.getMass() + " value: " + selectedPeak.getValue());
                            repaint();
                        }
                    }
                }
            }
        });
        this.addKeyListener(new KeyAdapter() {
            @Override
            public void keyReleased(KeyEvent e) {
                char keyChar = Character.toLowerCase(e.getKeyChar());
                if (keyChar == '+') {
                    scale *= 1.1;
                    updateDimension();
                }
                if (keyChar == '-') {
                    scale *= 0.9;
                    updateDimension();
                }
            }
        });
    }

    private TooltipCandidate getTooltip(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
        for (TooltipCandidate tooltipCandidate : tooltips) {
            if (tooltipCandidate.isValid(x, y / LINE_HEIGHT + 1, scale)) {
                return tooltipCandidate;
            }
        }
        return null;
    }

    public void setProtein(Protein protein) {
        this.proteinId = protein.getProteinId();
        this.protein = protein;
        proteinSpectrum = ShiftEngine.getSpectrum(protein.getSimplifiedAcids());
        if (scan != null) {
            initBestShift();
        } else {
            bestShiftB = 0;
            bestShiftY = 0;
        }
    }

    public Scan createReducedScan() {
        if (protein != null && scan != null) {
            List<Peak> reducedPeaks = new ArrayList<Peak>();
            for (Peak peak : scan.getPeaks()) {
                if (getPeakColor(peak) == Color.BLACK &&
                        getPeakColor(peak.getYPeak()) == Color.BLACK) {
                    reducedPeaks.add(peak);
                }
            }
            return new Scan(scan, reducedPeaks, proteinId);
        }
        return null;
    }

    public Scan getScan() {
        return scan;
    }

    private int totalTags;

    public boolean setScan(Scan scan) {
        if (scan != this.scan) {
            this.scan = scan;
            Analyzer analyzer = new Analyzer(conf);
            double precursorMassShift = PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan);
            peaks = scan.createSpectrumWithYPeaks(precursorMassShift);
            List<List<Peak>> components = analyzer.getComponents(peaks);
            this.components.clear();
            totalTags = 0;
            for (List<Peak> component : components) {
                List<Peak[]> tags = analyzer.getTags(component);
                totalTags += tags.size();
                this.components.add(tags);
            }

            if (proteinSpectrum != null) {
                initBestShift();
            }

            updateDimension();
            return true;
        }
        return false;
    }

    private void updateDimension() {
        double x = scan.getPrecursorMass();
        if (proteinSpectrum != null) {
            double xp = proteinSpectrum[proteinSpectrum.length - 1];
            if (xp > x) {
                x = xp;
            }
        }
        dimension = new Dimension((int) (x * scale) + 400, (components.size() + totalTags + 6) * LINE_HEIGHT);
        revalidate();
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
        int line = 1;

        if (scan == null)
            return;

        for (Peak peak : peaks) {
            double shift = peak.getPeakType() == PeakType.B ? bestShiftB : bestShiftY;
            drawPeak(g, peak, line, -shift);
        }
        double precursorMass = scan.getPrecursorMass();
        double end = precursorMass + bestShiftB;
        drawPeak(g, line, bestShiftB);
        drawPeak(g, line, end);
        g.drawString(df.format(precursorMass) + " + " + df.format(bestShiftB) + " " + df.format(bestShiftY), (int) (end * scale) + 3, 15);

        line++;

        if (protein != null) {
            String sequence = protein.getSimplifiedAcids();
            double pos = 0;
            for (int cur = 0; cur < sequence.length(); cur++) {
                pos += drawLetter(g, line, pos, Acid.getAcid(sequence.charAt(cur)));
                drawPeak(g, line, pos);
            }
            line++;

            int delta = 0;

            List<MassMatch> goodMatchesCopy = new ArrayList<MassMatch>();
            goodMatchesCopy.addAll(goodMatches);
            int[] stat = new int[proteinSpectrum.length];
            for (MassMatch match : goodMatches) {
                stat[match.getBegin()]++;
                stat[match.getEnd()]++;
            }


            while (goodMatchesCopy.size() > 0) {
                delta = drawBestCleavage(g, goodMatchesCopy, stat, line, delta);
                delta += 2;
            }

            if (bestDiff != 0) {
                g.setColor(Color.RED);
                for (MassMatch badMatch : badMatches) {
                    double mass = badMatch.getMass();
                    for (int i = 0; i < proteinSpectrum.length; i++) {
                        for (int j = i + 1; j < proteinSpectrum.length; j++) {
                            MassMatch nextMatch = new MassMatch(mass, i, j, bestDiff, proteinSpectrum);
                            if (nextMatch.getRelativeError() < PlaceStatistics.TEN_PPM) {
                                int y = line * LINE_HEIGHT + delta;
                                int x1 = (int) (proteinSpectrum[i] * scale);
                                int x2 = (int) (proteinSpectrum[j] * scale);
                                g.drawLine(x1, y, x2, y);
                                delta +=2;
                            }
                        }
                    }
                }
                g.setColor(Color.BLACK);
            }

            line += (delta + 2 * LINE_HEIGHT)/LINE_HEIGHT;

            for (Double shift : rangeShifts) {
                drawSpectrumAgainstProtein(g, line, precursorMass, shift);
                line++;
            }

        }

        for (int i = 0; i < components.size(); i++) {
            List<Peak[]> component = components.get(i);
            double min = scan.getPrecursorMass();
            for (Peak[] tag : component) {
                double v = tag[0].getValue();
                if (v < min) {
                    min = v;
                }
            }
            g.drawString("Component " + (i + 1), 3, line * LINE_HEIGHT);
            line++;
            for (Peak[] tag : component) {
                for (int j = 0; j < tag.length; j++) {
                    Peak peak = tag[j];
                    double peakValue = drawPeak(g, peak, line, min);
                    if (j + 1 < tag.length) {
                        double delta = tag[j + 1].getValue() - peakValue;
                        drawLetter(g, line, peakValue - min, Acid.getAcid(delta));
                    }
                }
                line++;
            }
        }
    }

    private int drawBestCleavage(Graphics g, List<MassMatch> goodMatches, int[] stat, int line, int delta) {
        int bestScore = -1;
        int bestPos = 0;
        for (int i = 0; i < stat.length; i++) {
             if (stat[i] > bestScore) {
                 bestPos = i;
                 bestScore = stat[i];
             }

        }

        for (Iterator<MassMatch> it = goodMatches.iterator(); it.hasNext(); ) {
            MassMatch match =  it.next();
            if (match.getBegin() == bestPos || match.getEnd() == bestPos) {
                double modAbs = Math.abs(match.getMod());
                g.setColor(Color.BLACK);
                if (modAbs == 0) {
                    g.setColor(Color.GREEN);
                }

                int y = line * LINE_HEIGHT + delta;
                int x1 = (int) (proteinSpectrum[match.getBegin()] * scale);
                int x2 = (int) (proteinSpectrum[match.getEnd()] * scale);
                int d = (x2 - x1) / 4;
                double mod = match.getMod();

                int x1mod = mod > 0 ? x1 : x1 + d + 2;
                int x2mod = mod > 0 ? x2 - d - 2 : x2;
                if (mod == 0) {
                    x1mod = x1;
                }
                g.drawLine(x1mod, y, x2mod, y);



                if (modAbs == 1) {
                    g.setColor(Color.ORANGE);
                }
                if (modAbs == Consts.WATER) {
                    g.setColor(Color.BLUE);
                }
                if (modAbs == Consts.AMMONIA) {
                    g.setColor(Color.MAGENTA);
                }


                if (mod < 0) {
                    g.drawLine(x1, y, x1 + d, y);
                } else if (mod > 0) {
                    g.drawLine(x2 - d, y, x2, y);
                }
                delta+=2;
                it.remove();
            }
        }
        stat[bestPos] = 0;
        return delta;
    }

    private void drawSpectrumAgainstProtein(Graphics g, int line, double precursorMass, double shift) {
        double end;
        g.setColor(Color.GREEN);
        for (Peak peak : peaks) {
            drawSharedPeak(g, peak, line, shift);
        }
        end = precursorMass + shift;
        updateColor(g, shift);
        drawPeak(g, line, shift);
        updateColor(g, end);
        drawPeak(g, line, end);
        g.setColor(Color.BLACK);
        g.drawString(df.format(shift), 20, LINE_HEIGHT * line);
    }

    private void updateColor(Graphics g, double value) {
        g.setColor(ShiftEngine.contains(proteinSpectrum, value) ? Color.GREEN : Color.BLACK);
    }

    private double drawPeak(Graphics g, Peak peak, int y, double min) {
        double peakValue = peak.getValue();
        double value = (peakValue - min);
        if (proteinSpectrum != null) {
            Color color = getPeakColor(peak);
            g.setColor(color);
        }
        drawPeak(g, y, value, peak);
        g.setColor(Color.BLACK);
        return peakValue;
    }

    private double drawSharedPeak(Graphics g, Peak peak, int y, double shift) {
        double value = peak.getValue() + shift;
        if (ShiftEngine.contains(proteinSpectrum, value)) {
            drawPeak(g, y, value, peak);
        }
        return value;
    }

    private Color getPeakColor(Peak peak) {
        Color color = Color.BLACK;
        double peakValue = peak.getValue();
        double shift = peak.getPeakType() == PeakType.B ? bestShiftB : bestShiftY;
        double v = peakValue + shift;
        if (ShiftEngine.contains(proteinSpectrum, v)) {
            color = Color.GREEN;
        } else if (ShiftEngine.contains(proteinSpectrum, v,
                v - 1, v + 1,
                v + Consts.WATER, v - Consts.WATER, v - Consts.AMMONIA, v + Consts.AMMONIA)) {
            color = Color.ORANGE;
        }
        return color;
    }

    private void initBestShift() {
        System.out.println("+ PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan); = " + +PrecursorMassShiftFinder.getPrecursorMassShift(conf, scan));
        System.out.print("Calculate best shift for B-ions...");
        bestShiftB = ShiftEngine.getBestShift(scan.getPeaks(), proteinSpectrum);
        System.out.println(Double.toString(bestShiftB));
        System.out.print("Calculate best shift for Y-ions...");
        bestShiftY = ShiftEngine.getBestShift(scan.getYPeaks(), proteinSpectrum);
        System.out.println(Double.toString(bestShiftY));
        rangeShifts = ShiftEngine.getBestShifts(scan.getPeaks(), scan.getPrecursorMass(), proteinSpectrum);
        List<MassMatch> matches = MassAlign.getMassAlign(proteinSpectrum, scan.getPeaks());
        goodMatches = new ArrayList<MassMatch>();
        badMatches = new ArrayList<MassMatch>();
        for (MassMatch match : matches) {
            if (match.getRelativeError() < PlaceStatistics.TEN_PPM) {
                goodMatches.add(match);
            } else {
                badMatches.add(match);
            }
        }

        bestDiff = ModificationDiscoverer.getBestDiff(scan.getId(), proteinId, proteinSpectrum, badMatches);
        String text = goodMatches.size() + " aligned masses out of " + scan.getPeaks().size() + " peaks.";
        if (bestDiff != 0) {
             text += " Best modification is " + bestDiff;
        }
        tagFinder.updateStatus(text);
    }

    private void drawPeak(Graphics g, int y, double value) {
        drawPeak(g, y, value, null);
    }

    private void drawPeak(Graphics g, int line, double value, Peak peak) {
        int v = (int) (value * scale);
        int y = line * LINE_HEIGHT;
        g.drawLine(v, y, v, y - LINE_HEIGHT + 3);
        double tooltipValue = peak != null ? peak.getValue() : value;
        tooltips.add(new TooltipCandidate(value, line, df.format(tooltipValue), peak));
        if (peak != null) {
            drawIon(g, line, v, peak.getPeakType());
            if (peak == selectedPeak) {
                Color prevColor = g.getColor();
                g.setColor(Color.BLUE);
                g.drawRect(v - 5, y - LINE_HEIGHT, 10, LINE_HEIGHT + 3);
                g.setColor(prevColor);
            }
        }
    }

    private void drawIon(Graphics g, int line, int v, PeakType peakType) {
        int x1 = peakType == PeakType.B ? v - 3 : v + 3;
        int y = line * LINE_HEIGHT;
        g.drawLine(x1, y - LINE_HEIGHT + 3, v, y - LINE_HEIGHT + 3);
        g.drawLine(x1, y, v, y);
    }

    private double drawLetter(Graphics g, int line, double pos, Acid acid) {
        double delta = acid.getMass();
        double x = pos + delta / 2 - 3;
        x *= scale;
        g.drawString(acid.name(), (int) (x), line * LINE_HEIGHT - 4);
        return delta;
    }
}
