package ru.spbau.bioinf.tagfinder.ui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.Map;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.KeyStroke;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.EValueAdapter;
import ru.spbau.bioinf.tagfinder.Protein;

import java.util.List;
import ru.spbau.bioinf.tagfinder.Scan;

public class TagFinder extends JFrame {

    private Configuration conf;
    private List<Protein> proteins;
    private Map<Integer,Integer> msAlignResults;
    private Map<Integer,Scan> scans;
    private JTabbedPane tabs;
    private List<JPanel> tabsList = new ArrayList<JPanel>();

    public TagFinder(String[] args) throws Exception {
        super("TagFinder");
        conf = new Configuration(args);
        proteins = conf.getProteins();
        msAlignResults = conf.getMSAlignResults();
        EValueAdapter.init(conf);

        tabs = new JTabbedPane();
        scans = conf.getScans();
        JPanel scanPanel = new ScanPanel(conf, scans, proteins, msAlignResults);
        addTab("Scan", scanPanel);
        JPanel proteinPanel = new ProteinPanel(proteins);
        addTab("Protein", proteinPanel);
        this.getContentPane().add(tabs);
        JMenuBar menubar = new JMenuBar();
        JMenu scanMenu = new JMenu("Scan");
        scanMenu.setMnemonic(KeyEvent.VK_S);
        JMenuItem evalue = new JMenuItem("E-Value");
        evalue.setMnemonic(KeyEvent.VK_E);
        evalue.setAccelerator(KeyStroke.getKeyStroke(
                KeyEvent.VK_E, ActionEvent.CTRL_MASK));


        evalue.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JPanel panel = tabsList.get(tabs.getSelectedIndex());
                if (panel instanceof ScanPanel) {
                    ((ScanPanel) panel).calculateEValue();
                }
            }
        });
        JMenuItem reduce = new JMenuItem("Reduce");
        reduce.setMnemonic(KeyEvent.VK_R);
        reduce.setAccelerator(KeyStroke.getKeyStroke(
                KeyEvent.VK_R, ActionEvent.CTRL_MASK));
        reduce.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JPanel panel = tabsList.get(tabs.getSelectedIndex());
                if (panel instanceof ScanPanel) {
                    Scan reducedScan = ((ScanPanel) panel).createReducedScan();
                    if (reducedScan != null) {
                        addScanTab(reducedScan);
                    }
                }
            }
        });

        scanMenu.add(evalue);
        scanMenu.add(reduce);

        menubar.add(scanMenu);

        setJMenuBar(menubar);
    }

    private void addTab(String name, JPanel scanPanel) {
        tabs.addTab(name, scanPanel);
        tabsList.add(scanPanel);
    }

    public void addScanTab(Scan scan) {
        ScanPanel scanPanel = new ScanPanel(conf, scans, proteins, msAlignResults);
        scanPanel.setScan(scan);
        addTab(scan.getName(), scanPanel);
    }

    public static void main(String[] args) throws Exception {
        TagFinder frame = new TagFinder(args);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(1000, 500);
        frame.setVisible(true);
    }

}
