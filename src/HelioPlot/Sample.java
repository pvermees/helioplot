package HelioPlot;

import Jama.Matrix;

public class Sample {
    
    /* Sample(double[] row)
     * constructs a sample from a row of data containing the following elements:
     * U, sU, Th, sTh, Sm, sSm, He, sHe, C
     */
    public Sample(double[] row){
        this.sample = row;
    }
    
    /* double[][] cov(boolean doSm)
     * returns the [nxn] covariance matrix
     */
    public Matrix cov(boolean doSm) throws Exception {
        int n = doSm ? 3 : 2;
        double[][] cov = new double[n][n];
        cov[0][0] = vV();
        cov[0][1] = covVW();
        cov[1][0] = cov[0][1];
        cov[1][1] = vW();
        if (n==3){
            cov[0][2] = cov[0][1];            
            cov[1][2] = cov[0][1];
            cov[2][0] = cov[0][1];
            cov[2][1] = cov[0][1];
            cov[2][2] = vX();
        }
        Matrix out = new Matrix(cov);
        return out;
    }
    
    public double V() throws Exception {
        return Math.log10(U()/He());
    }

    public double vV() throws Exception {
        return  ((sU()/U())*(sU()/U()) + (sHe()/He())*(sHe()/He())) / (l10*l10);
    }
    
    public double sV() throws Exception {
        return Math.sqrt(vV());
    }
        
    public double W() throws Exception {
        return Math.log10(Th()/He());
    }

    public double vW() throws Exception {
        return ((sTh()/Th())*(sTh()/Th()) + (sHe()/He())*(sHe()/He())) / (l10*l10);
    }

    public double sW() throws Exception {
        return Math.sqrt(vW());
    }
    
    public double X() throws Exception {
        return Math.log10(Sm()/He());
    }

    public double vX() throws Exception {
        return ((sSm()/Sm())*(sSm()/Sm()) + (sHe()/He())*(sHe()/He())) / (l10*l10);
    }

    public double sX() throws Exception {
        return Math.sqrt(vX());
    }

    public double covVW() throws Exception {
        return ((sHe()/He())*(sHe()/He())) / (l10*l10);
    }
    
    public double U() throws Exception {
        return sample[0];
    }

    public double sU() throws Exception {
        return sample[1];
    }

    public double Th() throws Exception {
        return sample[2];
    }
    
    public double sTh() throws Exception {
        return sample[3];
    }
    
    public double Sm() throws Exception {
        return sample[4];
    }

    public double sSm() throws Exception {
        return sample[5];
    }
    
    public double He() throws Exception {
        return sample[6];
    }
    
    public double sHe() throws Exception {
        return sample[7];
    }  
    
    public double C() throws Exception {
        return sample[8];
    }
    
    protected double[] sample;
    protected double l10 = Math.log(10);
}