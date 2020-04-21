/*
 * Based on Sun's ColorChooserDemo
 */

package HelioPlot;

import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.colorchooser.*;

/* ColorChooserDemo.java requires no other files. */
public class ColourChooser extends JPanel
                              implements ChangeListener {

    private ColourChooser(Coloured colouredobject) {
        super(new BorderLayout());
        this.colouredobject = colouredobject;

        // Set up the banner at the top of the window
        banner = new ColourBanner(colouredobject);
        
        // Add the radio buttons     
        minRadioButton = new javax.swing.JRadioButton();
        maxRadioButton = new javax.swing.JRadioButton();        
        minRadioButton.setSelected(true);
        minRadioButton.setText("min");
        maxRadioButton.setText("max");
        ButtonGroup minmax = new ButtonGroup();
        minmax.add(this.minRadioButton);
        minmax.add(this.maxRadioButton);   
        
        JPanel bannerPanel = new JPanel(new BorderLayout());
        bannerPanel.add(minRadioButton, BorderLayout.BEFORE_LINE_BEGINS);
        bannerPanel.add(banner, BorderLayout.CENTER);
        bannerPanel.add(maxRadioButton, BorderLayout.AFTER_LINE_ENDS);
        bannerPanel.setBorder(BorderFactory.createTitledBorder("Preview"));

        //Set up color chooser
        int[] minRGB = colouredobject.getMinRGB();
        Color initColor = new Color(minRGB[0],minRGB[1],minRGB[2]);
        tcc = new JColorChooser(initColor);
        tcc.setBorder(BorderFactory.createTitledBorder("Pick a Colour"));
        tcc.setPreviewPanel(new JPanel()); // remove preview pane
        this.removePanes();
        
        add(bannerPanel, BorderLayout.CENTER);
        add(tcc, BorderLayout.PAGE_END);       
    }
    
    public static ColourChooser createColourChooser(Coloured colouredobject) {
        ColourChooser instance = new ColourChooser(colouredobject);
        instance.tcc.getSelectionModel().addChangeListener(instance);
        return instance;
    }
    
    private void removePanes(){
        // Retrieve the current set of panels
        AbstractColorChooserPanel[] oldPanels = tcc.getChooserPanels();
        for (AbstractColorChooserPanel oldPanel : oldPanels) {
            String clsName = oldPanel.getClass().getName();
            if (clsName.equals("javax.swing.colorchooser.DefaultSwatchChooserPanel")) {
            } else if (clsName.equals("javax.swing.colorchooser.DefaultRGBChooserPanel")) {
            } else if (clsName.equals("javax.swing.colorchooser.DefaultHSBChooserPanel")) {
                // Remove hsb chooser if desired
                tcc.removeChooserPanel(oldPanel);
            }
        } 
    }

    public void stateChanged(ChangeEvent e) {
        Color newColor = tcc.getColor();
        this.refresh(newColor);
    }
    
    private void refresh(Color c){
        int r = c.getRed(), g = c.getGreen(), b = c.getBlue();
        if (this.minRadioButton.isSelected()){
            colouredobject.setMinRGB(r, g, b);
        } else if (this.maxRadioButton.isSelected()){
            colouredobject.setMaxRGB(r, g, b);
        }
        banner.repaint();
    }     

    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
     * @param colouredobject
     */
    public static void createAndShowGUI(ColouredObject colouredobject) {
        //Create and set up the window.
        final JFrame frame = new JFrame("Colour Chooser");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        //Create and set up the content pane.
        JComponent newContentPane = ColourChooser.createColourChooser(colouredobject);
        newContentPane.setOpaque(true); //content panes must be opaque
        frame.setContentPane(newContentPane);

        // make sure that the colour chooser is in focus until you close it
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
    
    Coloured colouredobject;
    protected JColorChooser tcc;
    protected ColourBanner banner;
    protected JRadioButton minRadioButton, maxRadioButton;
    
}