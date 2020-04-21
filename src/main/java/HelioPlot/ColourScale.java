package HelioPlot;

import java.awt.Color;
import java.awt.Graphics2D;
import java.text.DecimalFormat;

public class ColourScale extends ColouredObject {

    public ColourScale(HelioPlot plot){
        super();
        try {
            this.plot = plot;
            this.setMinMaxC();
            this.setMinRGB(255, 0, 255);
            this.setMaxRGB(0, 255, 0);
        } catch (Exception ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }

    final public void setMinMaxC() throws Exception {
        try {
            double[] minmaxc = plot.getMinMaxC();
                     this.minC = minmaxc[0];
                     this.maxC = minmaxc[1];
        } catch (Exception ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }

    @Override
    public void plot(){
        try {
            double[] minmaxc = plot.getMinMaxC();
            this.minC = minmaxc[0];
            this.maxC = minmaxc[1];
            this.initColours();
            this.plotBox();
            this.plotLabels();
        } catch (Exception ex){
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }
    
    private void plotLabels() throws Exception {
        DecimalFormat f = new DecimalFormat("#.##");
        float x1 = (float)(plot.wmap(plot.leftmargin)),
              x2 = x1 + (float)((1 - plot.leftmargin - plot.rightmargin)*
                   plot.getPanelWidth())-
                   (int)(plot.getGraphics2D().getFontMetrics().stringWidth(f.format(this.maxC))),
              y = (float) (plot.hmap(plot.bottommargin/15));
        plot.getGraphics2D().drawString(f.format(this.minC),x1,y);
        plot.getGraphics2D().drawString(f.format(this.maxC),x2,y);
        plot.getGraphics2D().drawString(plot.getColourLabel(),(x1+x2)/2,y);
    }
    
    private double getYpos(){
        return 2*plot.bottommargin/3;
    }
    
    private int getScaleHeight(){
        return (int)((1-plot.bottommargin)*plot.getPanelHeight()/15);
    }
    
    private void plotBox() throws Exception {
        Color c;
        Graphics2D g2 = plot.getGraphics2D();
        int x = plot.wmap(plot.leftmargin),
            x0 = x,
            y = plot.hmap(getYpos()),
            w = (int)((1 - plot.leftmargin - plot.rightmargin)*
                      plot.getPanelWidth())/this.numcolours,
            h = getScaleHeight();
        for (int i=0; i<this.numcolours; i++){
            c = new Color(colours[0][i],colours[1][i],colours[2][i]);
            g2.setColor(c);
            g2.fillRect(x, y, w, h);
            x += w;
            if (i+1<numcolours){
                w = (int)(((1 - plot.leftmargin - plot.rightmargin)*
                          plot.getPanelWidth()-(x-x0))/(this.numcolours-i-1));
            }
        }
        g2.setColor(Color.BLACK);
        x = plot.wmap(plot.leftmargin);
        y = plot.hmap(getYpos());
        w = (int)((1 - plot.leftmargin - plot.rightmargin)*plot.getPanelWidth());
        g2.drawRect(x,y,w,h);
    }
    
    private HelioPlot plot;

}