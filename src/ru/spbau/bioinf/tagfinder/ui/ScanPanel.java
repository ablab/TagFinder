package ru.spbau.bioinf.tagfinder.ui;


import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.Protein;
import ru.spbau.bioinf.tagfinder.Scan;

public class ScanPanel extends JPanel {

    private Map<Integer, Scan> scans;

    int scanId = 0;

    private final AtomicBoolean needUpdate = new AtomicBoolean(false);

    private final JLabel proteinName = new JLabel();

    private final JLabel proteinLabel = new JLabel("Protein ID: ");
    private final JTextField proteinIdInput = new JTextField();

    private final JLabel scanLabel = new JLabel("Scan ID: ");
    private final JLabel scanIdValueLabel = new JLabel();


    private final JLabel scanIdInputLabel = new JLabel("Enter new scan ID: ");
    private final JTextField scanIdInput = new JTextField();

    private ScanView scanView;
    private Map<Integer,Integer> msAlignResults;

    public ScanPanel(Configuration conf, Map<Integer, Scan> scans, List<Protein> proteins, Map<Integer,Integer> msAlignResults) {
        this.msAlignResults = msAlignResults;
        scanView = new ScanView(conf, proteins);
        this.scans = scans;
        scanIdInput.addKeyListener(new KeyAdapter() {
            @Override
            public void keyTyped(KeyEvent keyEvent) {
                checkNewScanId();
            }

            @Override
            public void keyPressed(KeyEvent keyEvent) {
                checkNewScanId();
            }

            @Override
            public void keyReleased(KeyEvent keyEvent) {
                checkNewScanId();
            }
        });
        scanIdInput.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                checkNewScanId();
            }
        });

        GridBagLayout gbl = new GridBagLayout();
        this.setLayout(gbl);
        GridBagConstraints gbc = new GridBagConstraints();

        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 0;
        add(scanLabel, gbc);
        gbc.gridx++;
        gbc.weightx = 1;
        gbc.anchor = GridBagConstraints.LINE_START;
        add(scanIdValueLabel, gbc);

        gbc.gridx++;
        gbc.weightx = 0;
        gbc.anchor = GridBagConstraints.LINE_END;
        add(proteinLabel, gbc);

        gbc.gridx++;
        gbc.weightx = 1;
        gbc.anchor = GridBagConstraints.LINE_START;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        proteinIdInput.setColumns(10);
        add(proteinIdInput, gbc);

        proteinIdInput.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                try {
                    int newProteinId = Integer.parseInt(proteinIdInput.getText());
                    needUpdate.compareAndSet(false, scanView.setProteinId(newProteinId));
                    update();
                } catch (NumberFormatException e1) {
                    //just wrong input
                }
            }
        });

        gbc.gridx++;
        gbc.weightx = 0;
        gbc.anchor = GridBagConstraints.LINE_END;
        add(scanIdInputLabel, gbc);

        gbc.gridx++;
        gbc.weightx = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.LINE_START;
        scanIdInput.setColumns(7);
        add(scanIdInput, gbc);


        gbc.gridx = 0;
        gbc.gridy++;
        gbc.gridwidth = 6;
        gbc.weighty = 0;
        gbc.anchor = GridBagConstraints.LINE_START;
        proteinName.setText("Protein name: ");
        add(proteinName, gbc);

        gbc.fill = GridBagConstraints.BOTH;
        gbc.gridy++;
        gbc.weighty = 1;

        JScrollPane scrollScanView = new JScrollPane(scanView);
        scrollScanView.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollScanView.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);

        add(scrollScanView, gbc);
        needUpdate.set(true);
        update();
    }

    private void checkNewScanId() {
        try {
            int newScanId = Integer.parseInt(scanIdInput.getText());
            if (newScanId != scanId) {
                scanId = newScanId;
                Scan scan = scans.get(scanId);
                if (scan != null) {
                    needUpdate.compareAndSet(false, scanView.setScan(scan));
                    update();
                }
            }
        } catch (NumberFormatException e) {
            //Nothing special - just text in number field;
        }
    }

    public void update() {
        if (needUpdate.getAndSet(false)) {
            Protein protein = scanView.getProtein();
            if (protein != null) {
                proteinName.setText("Protein " + protein.getProteinId() + " " + protein.getName());
            }

            if (msAlignResults.containsKey(scanId)) {
                int proteinId = msAlignResults.get(scanId);
                proteinLabel.setText("Protein ID (" + proteinId + ")");
            }

            scanIdValueLabel.setText(Integer.toString(scanId));
            repaint();
        }
    }
}
