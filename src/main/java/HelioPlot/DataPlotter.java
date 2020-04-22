package HelioPlot;

import Jama.Matrix;
import java.awt.Color;
import java.awt.Graphics2D;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;

abstract public class DataPlotter extends ColouredObject {

    public DataPlotter(HelioPlot plot){
        super();
        this.plot = plot;
    }

    protected void plotData() throws Exception {
        Data data = plot.getData();
        int ii = 0;
        boolean doSm = data.plotSm();
        for (Iterator<Sample> i = data.iterator(); i.hasNext(); ) {
            ii++;
            Sample sample = i.next();
            if (doSm){
                plotEllipse(new Ellipse(25,numsigma,sample),sample.Sm());
            } else {
                plotEllipse(new Ellipse(25,numsigma,sample),sample.C());
            }
            if (plot.plotNumber){
                plotNumber(ii,sample);
            }
        }
    }

    protected void plotLegend() throws Exception {
        DecimalFormat fm = new DecimalFormat("#.##");
        Graphics2D g2 = plot.getGraphics2D();
        String agestring;
        Data data = plot.getData();
        double[] ae;
        int numages = data.numAges();
        if (numages==1){
            Iterator<Sample> i = data.iterator();
            Sample sample = i.next();
            ae = Calculate.ageErr(sample);
        } else {
            ae = data.averageAgeErr(true);
        }
        agestring = "Arithmetic mean = " + Calculate.formatAgeLabel(ae[0],"#.00") +
                     ", s.e. = " + Calculate.formatAgeLabel(ae[1],"#.00");
        if (numages > 1){
            double mswdae = data.MSWD(ae[0],true);
            agestring += ", MSWD = " + fm.format(mswdae);
        }
        g2.drawString(agestring, plot.wmap(plot.leftmargin),
                      plot.hmap(1-plot.topmargin/4));
        if (numages==1){
            return;
        }
        double[] ge = data.averageAgeErr(false);
        double mswdge = data.MSWD(ge[0],false);
        String geomagestring = "Geometric mean = " + Calculate.formatAgeLabel(ge[0],"#.00") +
                    ",  s.e. = " + Calculate.formatAgeLabel(ge[1],"#.00") +
                    ", MSWD = " + fm.format(mswdge);
        g2.drawString(geomagestring, plot.wmap(plot.leftmargin),
                      plot.hmap(1-plot.topmargin/4)+2*plot.tlength);
        try {
            ae = plot.getData().centralAgeErr();
            double[] ci = data.CI(0.05);
            double mswdc = data.MSWD();
            agestring = "Central age = " + Calculate.formatAgeLabel(ae[0],"#.00");
            if (ae[1] != Data.NAN){
                     agestring += ", s.e. = " + Calculate.formatAgeLabel(ae[1],"#.00");
            }
            agestring += ", MSWD = " + fm.format(mswdc);
            String cistring = "95% C.I = [" + Calculate.formatAgeLabel(ci[0],"#.00") +
                              ", " + Calculate.formatAgeLabel(ci[1],"#.00") + "]";
            g2.drawString(agestring, plot.wmap(plot.leftmargin),
                     plot.hmap(1-plot.topmargin/4)+4*plot.tlength);
            g2.drawString(cistring, plot.wmap(plot.leftmargin),
                     plot.hmap(1-plot.topmargin/4)+6*plot.tlength);
        } catch (Exception ex){
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }

    // to plot the average logratio composition
    protected void plotAvgComp() throws Exception {
        Matrix[] MuCov = Newton.solveMuCov(plot.getData());
        plotEllipse(new Ellipse(25,numsigma,
                                   MuCov[0].get(0,0),Math.sqrt(MuCov[1].get(0,0)),
                                   MuCov[0].get(1,0),Math.sqrt(MuCov[1].get(1,1)),
                                   MuCov[1].get(0,1)),Data.NAN);
    }

    protected void plotEllipse(Ellipse ell, double c) throws Exception {
        Graphics2D g2 = plot.getGraphics2D();
        double[] xy;
        int[] x = new int[ell.size()],
              y = new int[ell.size()];
        for (int i=0; i<ell.size(); i++){
            xy = ellipse2xy(ell.getRow(i));
            x[i] = plot.wmap(xy[0]);
            y[i] = plot.hmap(xy[1]);
        }
        if (c == Data.NAN){
            g2.setColor(Color.white);
        } else {
            g2.setColor(plot.getColourScale().getColour(c));
        }
        g2.fillPolygon(x, y, ell.size());
        g2.setColor(Color.black);
        g2.drawPolygon(x, y, ell.size());
    }

    protected void plotNumber(int ii, Sample sample) throws Exception {
        Graphics2D g2 = plot.getGraphics2D();
        double[] xy = sample2xy(sample);
        g2.drawString(Integer.toString(ii), plot.wmap(xy[0]), plot.hmap(xy[1]));
    }
    
    abstract void autoScale() throws Exception;

    abstract double[] sample2xy(Sample sample) throws Exception ;

    abstract void plotContours() throws Exception ;

    abstract ArrayList<Double> getContourAges() throws Exception ;

    abstract double[] ellipse2xy(double[] row) throws Exception ;

    protected HelioPlot plot;
    protected final double NUMSIGMA = 2;
    protected double numsigma = NUMSIGMA;
}