package HelioPlot;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.text.DecimalFormat;
import javax.swing.*;

public class Settings extends JPanel {

    public Settings(Main parent){
        super();
        try {
            this.parent = parent;
            this.plot = parent.getHelioPlot();
            initMyComponents();
        } catch (Exception ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }

    private void initMyComponents(){
        initComponents();
        setValues();
    }
    
    public static void createAndShowGUI(Main parent) {
        //Create and set up the window.
        final JFrame frame = new JFrame("Settings");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        //Create and set up the content pane.
        JComponent newContentPane = new Settings(parent);
        newContentPane.setOpaque(true); //content panes must be opaque
        frame.setContentPane(newContentPane);

        // make sure that the frame is in focus until you close it
        frame.addWindowFocusListener(new WindowAdapter() {
            @Override
            public void windowLostFocus(WindowEvent evt) {
                frame.requestFocus();
            }
        });

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }    
    
    private void initComponents() {

        minV = new javax.swing.JTextField();
        jLabelMinV = new javax.swing.JLabel();
        maxV = new javax.swing.JTextField();
        jLabelMaxV = new javax.swing.JLabel();
        minW = new javax.swing.JTextField();
        jLabelMinW = new javax.swing.JLabel();
        maxW = new javax.swing.JTextField();
        jLabelMaxW = new javax.swing.JLabel();
        fU = new javax.swing.JTextField();
        jLabelfU = new javax.swing.JLabel();
        fTh = new javax.swing.JTextField();
        jLabelfTh = new javax.swing.JLabel();
        fHe = new javax.swing.JTextField();
        jLabelfHe = new javax.swing.JLabel();
        label = new javax.swing.JTextField();
        jLabelLabel = new javax.swing.JLabel();
        sigma = new javax.swing.JTextField();
        jLabelSigma = new javax.swing.JLabel();
        plotaverageRadioButton = new javax.swing.JRadioButton();
        numberRadioButton = new javax.swing.JRadioButton();
        OKbutton = new javax.swing.JButton();
        DefaultButton = new javax.swing.JButton();

        minV.setText("minV");

        jLabelMinV.setText("min log(U/He)");

        maxV.setText("maxV");

        jLabelMaxV.setText("max log(U/He)");

        minW.setText("minW");

        jLabelMinW.setText("min log(Th/He)");

        maxW.setText("maxW");

        jLabelMaxW.setText("max log(Th/He)");

        fU.setText("fU");

        jLabelfU.setText("x U");

        fTh.setText("fTh");

        jLabelfTh.setText("x Th");

        fHe.setText("fHe");

        jLabelfHe.setText("x He");

        label.setText("label");

        jLabelLabel.setText("label");

        sigma.setText("2");

        jLabelSigma.setText("σ errors");

        plotaverageRadioButton.setText("plot average composition");

        numberRadioButton.setText("plot sample numbers");

        OKbutton.setText("OK");
        OKbutton.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OKbuttonActionPerformed(evt);
            }
        });

        DefaultButton.setText("Defaults");
        DefaultButton.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                DefaultButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                            .addComponent(minV, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(minW, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(label, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabelMinW)
                            .addComponent(jLabelMinV)
                            .addComponent(jLabelLabel)))
                    .addComponent(numberRadioButton)
                    .addComponent(plotaverageRadioButton))
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                            .addComponent(sigma, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(maxW, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(maxV, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabelMaxV)
                            .addComponent(DefaultButton, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 79, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                            .addComponent(OKbutton, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(fHe, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(fTh, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(fU, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(20, 20, 20)
                                .addComponent(jLabelfU))
                            .addGroup(layout.createSequentialGroup()
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jLabelfHe)
                                    .addComponent(jLabelfTh, javax.swing.GroupLayout.Alignment.TRAILING)))))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(94, 94, 94)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabelSigma)
                            .addComponent(jLabelMaxW))))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        layout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {label, minV, minW});

        layout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {maxV, maxW, sigma});

        layout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {OKbutton, fHe, fTh, fU});

        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(minV, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabelMinV)
                    .addComponent(maxV, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabelMaxV)
                    .addComponent(fU, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabelfU))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(minW, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxW, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(fTh, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabelMinW)
                    .addComponent(jLabelMaxW)
                    .addComponent(jLabelfTh))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(label, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabelLabel)
                    .addComponent(sigma, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabelSigma)
                    .addComponent(fHe, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabelfHe))
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(plotaverageRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(numberRadioButton))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(27, 27, 27)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(DefaultButton)
                            .addComponent(OKbutton))))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
    }// </editor-fold>
    
    private void setValues() {
        try {
            TernaryPlot tplot = plot.getTernaryPlot();
            LogRatioPlot lplot = plot.getLogRatioPlot();
            fU.setText(fm.format(tplot.fU));
            fTh.setText(fm.format(tplot.fTh));
            fHe.setText(fm.format(tplot.fHe));
            minW.setText(fm.format(lplot.minW));
            minV.setText(fm.format(lplot.minV));
            maxV.setText(fm.format(lplot.maxV));
            maxW.setText(fm.format(lplot.maxW));
            sigma.setText(fm.format(lplot.numsigma));
            label.setText(plot.getLabel());
            numberRadioButton.setSelected(plot.plotNumber);
            plotaverageRadioButton.setSelected(plot.plotAverage);
        } catch (Exception ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }
    
    private void setDefaults() {
        try {
            plot.autoScale();
            TernaryPlot tplot = plot.getTernaryPlot();
            LogRatioPlot lplot = plot.getLogRatioPlot();            
            fU.setText(fm.format(tplot.fU));
            fTh.setText(fm.format(tplot.fTh));
            fHe.setText(fm.format(tplot.fHe));
            minW.setText(fm.format(lplot.minW));
            minV.setText(fm.format(lplot.minV));
            maxV.setText(fm.format(lplot.maxV));
            maxW.setText(fm.format(lplot.maxW));
            sigma.setText(fm.format(lplot.numsigma));
        } catch (Exception ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }

    private void DefaultButtonActionPerformed(java.awt.event.ActionEvent evt) {
        this.setDefaults();
    }    
    
    private void OKbuttonActionPerformed(java.awt.event.ActionEvent evt) {
        try{
            TernaryPlot tplot = plot.getTernaryPlot();
            LogRatioPlot lplot = plot.getLogRatioPlot();  
            lplot.minV = Double.valueOf(this.minV.getText());
            lplot.maxV = Double.valueOf(this.maxV.getText());
            lplot.minW = Double.valueOf(this.minW.getText());
            lplot.maxW = Double.valueOf(this.maxW.getText());
            tplot.fU = Double.valueOf(this.fU.getText());
            tplot.fTh = Double.valueOf(this.fTh.getText());
            tplot.fHe = Double.valueOf(this.fHe.getText());
            plot.setLabel(this.label.getText());
            plot.plotNumber = numberRadioButton.isSelected();
            plot.plotAverage = plotaverageRadioButton.isSelected();
            lplot.numsigma = Double.valueOf(this.sigma.getText());
            tplot.numsigma = lplot.numsigma;
            SwingUtilities.getWindowAncestor(this).dispose();
        } catch (Exception ex){
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
            this.setDefaults();
        }
    }
    
    // Variables declaration - do not modify
    private javax.swing.JButton DefaultButton;
    private javax.swing.JButton OKbutton;
    private javax.swing.JTextField fHe;
    private javax.swing.JTextField fTh;
    private javax.swing.JTextField fU;
    private javax.swing.JLabel jLabelLabel;
    private javax.swing.JLabel jLabelMaxV;
    private javax.swing.JLabel jLabelMaxW;
    private javax.swing.JLabel jLabelMinV;
    private javax.swing.JLabel jLabelMinW;
    private javax.swing.JLabel jLabelSigma;
    private javax.swing.JLabel jLabelfHe;
    private javax.swing.JLabel jLabelfTh;
    private javax.swing.JLabel jLabelfU;
    private javax.swing.JTextField label;
    private javax.swing.JTextField maxV;
    private javax.swing.JTextField maxW;
    private javax.swing.JTextField minV;
    private javax.swing.JTextField minW;
    private javax.swing.JRadioButton numberRadioButton;
    private javax.swing.JRadioButton plotaverageRadioButton;
    private javax.swing.JTextField sigma;
    // other private variables
    protected Main parent;
    protected HelioPlot plot;
    protected DecimalFormat fm = new DecimalFormat("#.###");
    // End of variables declaration

}