package ru.spbau.bioinf.tagfinder.ui;

import java.util.Map;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
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

    public TagFinder(String[] args) throws Exception {
        super("TagFinder");
        conf = new Configuration(args);
        proteins = conf.getProteins();
        msAlignResults = conf.getMSAlignResults();
        tabs = new JTabbedPane();
        scans = conf.getScans();
        JPanel scanPanel = new ScanPanel(conf, scans, proteins, msAlignResults, this);
        tabs.addTab("Scan", scanPanel);
        JPanel proteinPanel = new ProteinPanel(proteins);
        tabs.addTab("Protein", proteinPanel);
        this.getContentPane().add(tabs);

        EValueAdapter.init(conf);
    }

    public void addScanTab(Scan scan) {
        ScanPanel scanPanel = new ScanPanel(conf, scans, proteins, msAlignResults, this);
        scanPanel.setScan(scan);
        tabs.addTab(scan.getName(), scanPanel);
    }

    public static void main(String[] args) throws Exception {
        TagFinder frame = new TagFinder(args);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(1000, 500);
        frame.setVisible(true);
    }

}
