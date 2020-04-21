package HelioPlot;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import javax.swing.JPanel;

public class ColourBanner extends JPanel {
    
    public ColourBanner(Coloured colouredobject){
        this.colouredobject = colouredobject;
        this.setPreferredSize(new Dimension(100,65));    
    }
    
    @Override
    public void paintComponent(Graphics g){
        try {
            super.paintComponent(g);
            this.fillTheBanner(g);
        } catch (Exception ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }
    
    private void fillTheBanner(Graphics g) throws Exception {
        Color c;
        int[][] colours = colouredobject.getColours();
        int x = 0,
            numcolours = colouredobject.getNumColours(),
            w = (int)(this.getWidth()/numcolours),
            h = (int)(this.getHeight());
        for (int i=0; i<numcolours; i++){
            c = new Color(colours[0][i],colours[1][i],colours[2][i]);
            g.setColor(c);
            g.fillRect(x, 0, w, h);
            x += w;
            if (i+1<numcolours){
                w = (int)((this.getWidth()-x)/(numcolours-i-1));
            }
        }
        this.getGraphics().setColor(Color.BLACK);
        w = (int)(this.getWidth());
        this.getGraphics().drawRect(0,0,w,h);
    }

    Coloured colouredobject;
}
