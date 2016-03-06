package HelioPlot;

import Jama.Matrix;
import java.awt.Color;
import java.awt.Graphics2D;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;

public class TernaryPlot extends DataPlotter {

    public TernaryPlot(HelioPlot plot) {
        super(plot);
        try {
            this.setMinRGB(255, 255, 255);
            this.setMaxRGB(255, 0, 0);
        } catch (Exception ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }
    
    @Override
    public void plot(){
        try {
            plotContours();
            plotTriangle();
            plotData();
            plotLegend();
            if (plot.plotAverage){plotAvgComp();}
        } catch (Exception ex){
            if (Main.DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }
    
    @Override
    protected ArrayList<Double> getContourAges() throws Exception {
        double tmin, tmax, Sm = plot.getData().getSmBar();
        double[] pmin = this.zoomOut(0.95,0,0.05),
                 pmax = this.zoomOut(0.05,0,0.95);
        tmin = Calculate.age((1-Sm)*pmin[0], (1-Sm)*pmin[1], Sm, (1-Sm)*pmin[2]);
        tmax = Calculate.age((1-Sm)*pmax[0], (1-Sm)*pmax[1], Sm, (1-Sm)*pmax[2]);
        return getAgeContours(tmin,tmax,numcontours);
    }

    private void plotTriangle() throws Exception {
        Graphics2D g2 = plot.getGraphics2D();
        double[] XY1 = {0,0}, XY2 = {0.5,1}, XY3 = {1,0},
                 xy1 = plot.XY2xy(XY1), xy2 = plot.XY2xy(XY2), xy3 = plot.XY2xy(XY3);
        int[] xPoints = {plot.wmap(xy1[0]),plot.wmap(xy2[0]),plot.wmap(xy3[0])},
              yPoints = {plot.hmap(xy1[1]),plot.hmap(xy2[1]),plot.hmap(xy3[1])},
              maskx = {plot.wmap(0),xPoints[0],xPoints[1],plot.wmap(1),plot.wmap(0)},
              masky = {plot.hmap(0),yPoints[0],yPoints[1],plot.hmap(1),plot.hmap(1)};
        int dx, dy = g2.getFontMetrics().getHeight();
        g2.setColor(Color.WHITE);
        g2.fillPolygon(maskx, masky, 5);
        g2.setColor(Color.BLACK);
        g2.drawPolygon(xPoints, yPoints, 3);
        String[] labels = getLabels();
        g2.drawString(labels[0],xPoints[0], yPoints[0]+dy);
        dx = g2.getFontMetrics().stringWidth(labels[2])/2;
        g2.drawString(labels[2],xPoints[1]-dx, yPoints[1]-dy/2);
        dx = g2.getFontMetrics().stringWidth(labels[1]);
        g2.drawString(labels[1],xPoints[2]-dx, yPoints[2]+dy);        
    }
    
    private String[] getLabels() throws Exception {
        String[] labels = new String[3];
        DecimalFormat f = new DecimalFormat("0.##E0");
        labels[0] = (fU==1) ? "U" : String.format("%.2f", fU) + "xU";
        labels[1] = (fTh==1) ? "Th" : String.format("%.2f", fTh) + "xTh";
        labels[2] = (fHe==1) ? "He" : String.format("%.2f", fHe) + "xHe";
        return labels;
    }
    
    @Override
    double[] ellipse2xy(double[] row) throws Exception {
        return this.p2xy(row);
    }
    
    @Override
    public void plotContours() throws Exception {
        ArrayList<Double> t = getContourAges();
        this.minC = t.get(0);
        this.maxC = t.get(t.size()-1);
        this.maxC += (maxC-minC)/numcontours; // to ensure a different background
        Graphics2D g2 = plot.getGraphics2D();
        // fill triangle with maximum colour
        double[] XY1 = {0,0}, XY2 = {0.5,1}, XY3 = {1,0},
                 xy1 = plot.XY2xy(XY1), xy2 = plot.XY2xy(XY2), xy3 = plot.XY2xy(XY3);
        int[] xPoints = {plot.wmap(xy1[0]),plot.wmap(xy2[0]),plot.wmap(xy3[0])},
              yPoints = {plot.hmap(xy1[1]),plot.hmap(xy2[1]),plot.hmap(xy3[1])};
        g2.setColor(getColour(maxC));
        g2.fillPolygon(xPoints, yPoints, 3);
        // now get started with the contours
        double[] xy;
        int numsegments = 5;
        int[] x = new int[numsegments+2], y = new int[numsegments+2];
        Data data = plot.getData();
        double U, Th, Sm = data.getSmBar(); // if hasSm is false, then Sm is 0
        Collections.sort(t);
        // build segments
        for (int i=t.size()-1; i>=0; i--){
            // build the top part of the contour
            for (int j=0; j<numsegments; j++){
                Th = (j/(numsegments-1))*(1-Sm);
                U = Calculate.tThSm2U(t.get(i),Th,Sm);
                if (U < 0){
                    U = 0;
                    Th = Calculate.tUSm2Th(t.get(i),U,Sm);
                }
                if (1-U-Th-Sm < 0){
                    U = Calculate.tHeSm2U(t.get(i), 0, Sm);
                    Th = Calculate.tUSm2Th(t.get(i), U, Sm);
                }
                xy = UThSmt2xy(U, Th, Sm, t.get(i));
                x[j] = plot.wmap(xy[0]);
                y[j] = plot.hmap(xy[1]);
                if (xy[1]<plot.bottommargin) {
                    if (y[j]!=y[j-1]) x[j] = x[j-1]+(x[j]-x[j-1])*(plot.hmap(plot.bottommargin)-y[j-1])/(y[j]-y[j-1]);
                    y[j] = plot.hmap(plot.bottommargin);
                }
            }
            // close the contour so you can fill it with a colour
            x[numsegments] = plot.wmap(1-plot.rightmargin);
            y[numsegments] = plot.hmap(plot.bottommargin);
            x[numsegments+1] = plot.wmap(plot.leftmargin);
            y[numsegments+1] = plot.hmap(plot.bottommargin);
            // fill the contour
            g2.setColor(getColour(t.get(i)));
            g2.fillPolygon(x, y, numsegments+2);
            g2.setColor(Color.black);
            g2.drawPolygon(x, y, numsegments+2);
            // add label
            g2.drawString(Calculate.formatAgeLabel(t.get(i)),
                          x[numsegments-1], y[numsegments-1]);
        }
    }

    /* input: U, Th, Sm, t
     * output: XY = x, y coordinates for plotting on the ternary diagram
     * Sm is used for calculating t, but not for the ternary normalisation
     */
    private double[] UThSmt2xy(double U, double Th, double Sm, double t) throws Exception {
        double He;
        double[] p;
        He = Calculate.He(U, Th, Sm, t);
        p = Calculate.normalise(U,Th,He);
        return p2xy(p);
    }
    
    /* 
     * input: p = array of {p1,p2,p3} compositions
     * output: xy = array of x and y coordinates for plotting
     *              between 0 and 1, spanning the entire jPanel
     */
    private double[] p2xy(double[] p) throws Exception {
        double[] P = this.zoomIn(p);
        P = Calculate.normalise(P[0], P[1], P[2]);
        double[] XY = {0.5*(1-P[0]+P[1]),P[2]};
        double[] xy = plot.XY2xy(XY);
        return xy;
    } 

    @Override
    protected double[] sample2xy(Sample sample) throws Exception {
        double[] p = Calculate.normalise(sample.U(), sample.Th(), sample.He());
        return this.p2xy(p);
    }
    
    /* 
     * input: XY = array of X and Y coordinates between 0 and 1
     * output: p = array of {p1,p2,p3} compositions
     */
    public double[] XY2p(double[] XY) throws Exception {
        double he = XY[1],
               th = XY[0] - 0.5*XY[1],
               u = 1-th-he;
        double[] p = {u,th,he};
        //double[] p = Calculate.normalise(u, th, he);
        return this.zoomOut(p);
    }
    
    protected double[] zoomIn(double U, double Th, double He) throws Exception {
        double[] out = {U*fU,Th*fTh,He*fHe};
        return out;
    }    
    
    protected double[] zoomIn(double[] UThHe) throws Exception {
        return zoomIn(UThHe[0],UThHe[1],UThHe[2]);
    }
    
    protected double[] zoomOut(double u, double th, double he) throws Exception {
        double[] out = {u/fU,th/fTh,he/fHe};
        return out;
    } 
    
    protected double[] zoomOut(double[] p) throws Exception {
        return zoomOut(p[0],p[1],p[2]);
    }
    
    @Override
    public void autoScale() throws Exception {
        Data data = plot.getData();
        Matrix mle = data.mle2();
        double[] UThHe = Calculate.vw2UThHe(mle.getRowPackedCopy());
        double[] p = Calculate.normalise(UThHe[0], UThHe[1], UThHe[2]);
        fU = 1/p[0];
        fTh = 1/p[1];
        fHe = 1/p[2];
    }
        
    protected int numcontours = 5;
    protected final double FU = 1, FTH = 1, FHE = 100;
    protected double fU = FU, fTh = FTH, fHe = FHE;
        
}