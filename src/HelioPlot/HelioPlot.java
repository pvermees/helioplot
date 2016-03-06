package HelioPlot;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class HelioPlot extends JFrame {

    public HelioPlot(Main caller){
        super("Output");
        this.parent = caller;
        fixProportions();
        this.colourscale = new ColourScale(this);
        this.ternaryplot = new TernaryPlot(this);
        this.logratioplot = new LogRatioPlot(this);
        plotpanel = new JPanel(){
            @Override public void paintComponent(Graphics g){
                g2 = (Graphics2D) g;
                double scale = getScale();
                Font myfont = g2.getFont().deriveFont(AffineTransform.getScaleInstance(scale, scale));
                g2.setFont(myfont);
                tlength = (int)(0.5*g2.getFontMetrics().charWidth('W')/.75);
                super.paintComponent(g2);
                fixProportions();
                if (!parent.noColour()){
                    colourscale.plot();
                }
                if (parent.ternary()){
                    ternaryplot.plot();
                } else {
                    logratioplot.plot();
                }
            }
        };
        plotpanel.setBackground(Color.white);
        this.add(plotpanel);
        this.pack();
    }

    public int getPanelWidth(){
        return plotpanel.getWidth();
    }
    
    public int getPanelHeight(){
        return plotpanel.getHeight();
    }

    private double getScale(){
        double width = (double)this.getPanelWidth(),
               height = (double)this.getPanelHeight(),
               scalex = width/SIZE,
               scaley = height/SIZE;
        if (scalex<scaley){
            return scalex;
        } else {
            return scaley;
        }
    }

    public double[] getMinMaxC() throws Exception {
        Data data = parent.getData();
        ArrayList<Double> C = data.getC();
        double[] xy = new double[2];
        double min = Data.NAN, max = Data.NAN;
        int i = 0;
        if (data.plotSm()) {
            C = data.getSm();
        }
        // initialise
        for (; i<C.size(); i++){
            if (C.get(i)!=Data.NAN){
                min = C.get(i);
                max = C.get(i);
                break;
            }
        }
        // find min and max
        for (; i<C.size(); i++){
            if (C.get(i)!=Data.NAN && C.get(i)<min){min = C.get(i);}
            if (C.get(i)>max){max = C.get(i);}
        }
        xy[0] = min;
        xy[1] = max;
        return xy;
    }

    private void fixProportions(){
        try {
            bottommargin = BOTTOMMARGIN;
            rightmargin = RIGHTMARGIN;
            if ((1-bottommargin-topmargin)*getHeight()>
                Math.sqrt(0.75)*(1-leftmargin-rightmargin)*getPanelWidth()){
                bottommargin = 1 - topmargin -
                        ((1-leftmargin-rightmargin)*getPanelWidth()*Math.sqrt(0.75))/
                        getPanelHeight();
            } else {
                rightmargin = 1 - leftmargin -
                        (1-bottommargin-topmargin)*getPanelHeight()/
                        (Math.sqrt(0.75)*getPanelWidth());
            }
        } catch (Exception ex){
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
            bottommargin = BOTTOMMARGIN;
            rightmargin = RIGHTMARGIN;
        }
    }

    public int wmap(double val) throws Exception {
        return (int) (val*getPanelWidth());
    }

    public int hmap(double val) throws Exception {
        return (int) ((1-val)*getPanelHeight());
    }

    /*
     * This function pads the coordinates with margins for plotting
     * input: XY = plotContours coordinates (between 0 and 1)
     * output: xy = canvas coordinates (between 0 and 1)
     */
    public double[] XY2xy(double[] XY) throws Exception {
        double[] xy = new double[2];
        xy[0] = leftmargin + (1-leftmargin-rightmargin)*XY[0];
        xy[1] = bottommargin + (1-bottommargin-topmargin)*XY[1];
        return xy;
    }

    public void autoScale() throws Exception {
        ternaryplot.autoScale();
        logratioplot.autoScale();
    }

    public Graphics2D getGraphics2D(){
        return g2;
    }

    public ColourScale getColourScale(){
        return colourscale;
    }

    public TernaryPlot getTernaryPlot(){
        return ternaryplot;
    }

    public LogRatioPlot getLogRatioPlot(){
        return logratioplot;
    }

    public Data getData(){
        return parent.getData();
    }

    public void setLabel(String label){
        this.label = label;
    }

    public String getLabel(){
        return this.label;
    }
    
    public String getColourLabel() throws Exception{
        if (parent.data.plotSm()){
            return ("[Sm]");
        } else {
            return this.getLabel();
        }
    }

    @Override
    public Main getParent(){
        return parent;
    }

    JPanel plotpanel;
    private Graphics2D g2;
    private final double RIGHTMARGIN = 0.1, BOTTOMMARGIN = 0.15;
    double leftmargin = 0.1, rightmargin,
           bottommargin, topmargin = 0.2;
    protected int tlength = 6;
    private TernaryPlot ternaryplot;
    private LogRatioPlot logratioplot;
    private ColourScale colourscale;
    private Main parent;
    private String label = "[C]";
    protected boolean plotNumber = false, plotAverage = false;
    static final int SIZE = 500;
}