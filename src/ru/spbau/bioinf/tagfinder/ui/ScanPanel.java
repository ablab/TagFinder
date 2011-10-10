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
import javax.swing.JTextArea;
import javax.swing.JTextField;
import ru.spbau.bioinf.tagfinder.Protein;
import ru.spbau.bioinf.tagfinder.Scan;

public class ScanPanel extends JPanel {

    private List<Protein> proteins;
    private Map<Integer, Scan> scans;

    private int proteinId = 0;
    int scanId = 0;

    private final AtomicBoolean needUpdate = new AtomicBoolean(false);

    private final JLabel proteinName = new JLabel();
    private final JTextArea sequence = new JTextArea();

    private final JLabel proteinLabel = new JLabel("Protein ID: ");
    private final JLabel proteinIdValueLabel = new JLabel();

    private final JLabel scanLabel = new JLabel("Scan ID: ");
    private final JLabel scanIdValueLabel = new JLabel();

    private final JLabel tagTextLabel = new JLabel("Tag : ");
    private final JTextField tagInput = new JTextField();

    private final JLabel scanIdInputLabel = new JLabel("Enter new scan ID: ");
    private final JTextField scanIdInput = new JTextField();

    private ScanView scanView = new ScanView();

    public ScanPanel(Map<Integer, Scan> scans) {
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
        gbc.anchor = GridBagConstraints.LINE_END;
        gbc.weightx = 0;
        add(tagTextLabel, gbc);

        gbc.weightx = 1;
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.LINE_END;
        add(scanIdInputLabel, gbc);
        gbc.gridx++;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.LINE_START;
        scanIdInput.setColumns(7);
        add(scanIdInput, gbc);
        gbc.gridx = 0;
        gbc.gridy++;
        gbc.gridwidth = 5;
        gbc.weighty = 0;
        gbc.anchor = GridBagConstraints.LINE_START;
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
                needUpdate.set(true);
                update();
            }
        } catch (NumberFormatException e) {
            //Nothing special - just text in number field;
        }
    }

    public void update() {
        if (needUpdate.getAndSet(false)) {
            Scan scan = scans.get(scanId);
            if (scan != null) {
                scanView.setScan(scan);
                scanView.repaint();
                scanIdValueLabel.setText(Integer.toString(scanId));
            }
        }
    }
}
