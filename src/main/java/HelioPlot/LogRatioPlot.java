package HelioPlot;

import java.awt.Color;
import java.awt.Graphics2D;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

public class LogRatioPlot extends DataPlotter {

    public LogRatioPlot(HelioPlot plot){
        super(plot);
        this.setMinRGB(255, 255, 255);
        this.setMaxRGB(255, 0, 0);
    }

    @Override
    public void plot(){
        try {
            plotContours();
            plotAxes();
            plotData();
            plotLegend();
            if (plot.plotAverage){plotAvgComp();}
        } catch (Exception ex){
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }
    
    @Override
    protected ArrayList<Double> getContourAges() throws Exception {
        double[] vwmin = {this.minV,this.minW}, vwmax = {this.maxV,this.maxW};
        double[] xyzmin = Calculate.vw2UThHe(vwmin), xyzmax = Calculate.vw2UThHe(vwmax);
        double tmin = Calculate.age(xyzmax[0], xyzmax[1], 0, xyzmax[2]);
        double tmax = Calculate.age(xyzmin[0], xyzmin[1], 0, xyzmin[2]);
        return getAgeContours(tmin,tmax,numcontours);
    }

    private void plotAxes() throws Exception {
        Graphics2D g2 = plot.getGraphics2D();
        g2.drawLine(plot.wmap(plot.leftmargin), plot.hmap(plot.bottommargin), 
                    plot.wmap(1-plot.rightmargin), plot.hmap(plot.bottommargin));
        g2.drawLine(plot.wmap(plot.leftmargin), plot.hmap(plot.bottommargin), 
                    plot.wmap(plot.leftmargin), plot.hmap(1-plot.topmargin));
        g2.drawString("log(U/He)", plot.wmap(1-plot.rightmargin)+plot.tlength,
                                    plot.hmap(plot.bottommargin));
        g2.drawString("log(Th/He)", plot.wmap(plot.leftmargin)-
                                    (int)(0.5*g2.getFontMetrics().stringWidth("log(Th/He)")),
                                    plot.hmap(1-plot.topmargin)-2*plot.tlength);
        plotTicks(10);
    }

    private void plotTicks(int numticks) throws Exception {
        double vRange = this.maxV - this.minV;
        double wRange = this.maxW - this.minW;
        double dv, dw;
        DecimalFormat f = new DecimalFormat("#.##");
        if (vRange < 1 && wRange < 1) {
            dv = vRange/numticks;
            dw = wRange/numticks;
        } else {
            dv = 1;
            dw = 1;
        }
        Graphics2D g2 = plot.getGraphics2D();
        double[] xy, vw = {this.minV, this.minW};
        int i=0;
        // V
        while (true){
            vw[0] = this.minV + i*dv;
            if (vw[0]>this.maxV){break;}
            vw[1] = this.minW;
            xy = this.VW2xy(vw);
            g2.drawLine(plot.wmap(xy[0]), plot.hmap(xy[1]), 
                        plot.wmap(xy[0]), plot.hmap(xy[1])+plot.tlength);
            g2.drawString(f.format(vw[0]), plot.wmap(xy[0]), 
                          plot.hmap(xy[1])+ 
                          (g2.getFontMetrics().getHeight()));
            i++;
        }
        i=0;
        // W
        while (true){
            vw[0] = this.minV;
            vw[1] = this.minW + i*dw;
            if (vw[1]>this.maxW){break;}
            xy = this.VW2xy(vw);
            g2.drawLine(plot.wmap(xy[0]), plot.hmap(xy[1]), 
                        plot.wmap(xy[0]) - plot.tlength, plot.hmap(xy[1]));
            g2.drawString(f.format(vw[1]), 
                        plot.wmap(xy[0]) - plot.tlength -
                        (g2.getFontMetrics().stringWidth(f.format(vw[1]))),
                        plot.hmap(xy[1]));             
            i++;
        }
    }

    @Override
    double[] ellipse2xy(double[] row) throws Exception {
        return this.VW2xy(Calculate.UThHe2vw(row));
    }
    
    @Override
    public void plotContours() throws Exception {
        ArrayList<Double> t = getContourAges();        
        this.minC = t.get(0);
        this.maxC = t.get(t.size()-1);        
        Graphics2D g2 = plot.getGraphics2D();
        double[] xy = new double[2];
        double[] U_He = new double[numsegments];
        Data data = plot.getData();
        boolean doSm = data.hasSm();
        double foo, bar, U, Th, Sm = data.getSmBar(),
               SmHe = doSm ? Math.pow(10, data.mle3().get(2,0)) : 0, // average Sm/He ratio;
               He = doSm ? Sm/SmHe : 1;
        // build the segments in exponentially increasing U/He-steps
        for (int i=0; i<numsegments; i++) {
            foo = minV + i*(maxV-minV)/numsegments;
            bar = Math.pow(10,foo);
            U_He[i] = bar;
        }
        int[] x = new int[numsegments+2], y = new int[numsegments+2];
        // first two segments are always the same
        //x[0] = plot.wmap(1-plot.rightmargin);
        y[0] = plot.hmap(plot.bottommargin);
        x[1] = plot.wmap(plot.leftmargin);
        y[1] = plot.hmap(plot.bottommargin);
        Collections.sort(t);
        // loop through all the contours
        for (int i=0; i<t.size(); i++){
            // calculate the segment coordinates
            int labelx = 0, labely = plot.hmap(0);
            for (int j=0; j<numsegments; j++){
                U = He*U_He[j];
                xy = this.USmHet2xy(U, Sm, He, t.get(i)); 
                x[j+2] = plot.wmap(xy[0]);
                y[j+2] = plot.hmap(xy[1]);
                if (y[j+2]<=labely) { 
                    labely = y[j+2];
                    if (x[j+2]>labelx){
                        labelx = x[j+2];
                    }
                }
            }
            x[0] = plot.wmap(xy[0]);
            //y[0] = plot.hmap(plot.bottommargin);            
            // plot the contours as polygons
            g2.setColor(this.getColour(t.get(i)));
            g2.fillPolygon(x, y, numsegments+2);
            g2.setColor(Color.black);
            g2.drawPolygon(x, y, numsegments+2);
            // add label
            if (y[2]>=plot.hmap(plot.bottommargin)){
                // don't plot a label
            } else {
                g2.drawString(Calculate.formatAgeLabel(t.get(i)), labelx, labely);
            }
        }
    }

    /* input: U, Sm, He, t
     * output: xy = x, y coordinates for plotting
     */
    private double[] USmHet2xy(double U, double Sm, double He, double t) throws Exception {
        double[] xy, VW = new double[2];
        double W, Th = Calculate.USmHet2Th(U,Sm,He,t);
        VW[0] = Math.log10(U/He);
        if (Th<=0){
            VW[1] = minW;
        } else {
            W = Math.log10(Th/He);
            if (W>minW && W<maxW){
                VW[1] = W;
            } else if (W<minW) {
                VW[1] = minW;
            } else {
                VW[1] = maxW;
            }
        }
        xy = VW2xy(VW);
        return xy;
    }
    
    /* input: Th, Sm, He, t
     * output: xy = x, y coordinates for plotting
     */
    private double[] ThSmHet2xy(double Th, double Sm, double He, double t) throws Exception {
        double[] xy, VW = new double[2];
        double V, U = Calculate.ThSmHet2U(Th,Sm,He,t);
        VW[1] = Math.log10(Th/He);
        if (U<=0){
            VW[0] = minV;
        } else {
            V = Math.log10(U/He);
            if (V>minV && V<maxV){
                VW[0] = V;
            } else if (V<minV) {
                VW[0] = minV;
            } else {
                VW[0] = maxV;
            }
        }
        xy = VW2xy(VW);
        return xy;
    }
    
    /* input: VW = plotContours coordinates (between minV and maxV)
     * output: xy = canvas coordinates (between 0 and 1)
     */
    private double[] VW2xy(double[] VW) throws Exception {
        double[] XY = {(VW[0]-minV)/(maxV-minV), (VW[1]-minW)/(maxW-minW)};
        double[] xy = plot.XY2xy(XY);
        return xy;
    }
    
    @Override
    protected double[] sample2xy(Sample sample) throws Exception {
        double[] VW = Calculate.UThHe2vw(sample.U(), sample.Th(), sample.He());
        return this.VW2xy(VW);
    }
    
    @Override
    public void autoScale() throws Exception {
        double n = 3*numsigma;
        Data data = plot.getData();
        minV = MAXV;
        maxV = MINV;
        minW = MAXW;
        maxW = MINW;
        for (Iterator<Sample> i = data.iterator(); i.hasNext(); ) {
            Sample sample = i.next();
            if (sample.V()-n*sample.sV() < minV){
                minV = sample.V()-n*sample.sV();
            }
            if (sample.V()+n*sample.sV() > maxV){
                maxV = sample.V()+n*sample.sV();
            }
            if (sample.W()-n*sample.sW() < minW){
                minW = sample.W()-n*sample.sW();
            }
            if (sample.W()+n*sample.sW() > maxW){
                maxW = sample.W()+n*sample.sW();
            }
        }
    }    
    
    protected int numcontours = 8, numsegments = 98;
    protected final double FU = 1, FTH = 1, MINV = -1, MAXV = 6, MINW = -1, MAXW = 6;
    protected double fU = FU, fTh = FTH, minV = MINV, maxV = MAXV, minW = MINW, maxW = MAXW;    
    
}