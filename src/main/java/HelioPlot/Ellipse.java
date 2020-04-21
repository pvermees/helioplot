package HelioPlot;

public class Ellipse {

    /*
     * Ellipse(int n, double[] row) with 
     * n = resolution of the ellipse (default: 10)
     * row = 8-element vector of doubles (U, sU, Th, sTh, Sm, sSm, He, sHe)
     */
    public Ellipse(int n, double numsigma, Sample sample) throws Exception {
        double rho = sample.U()*sample.Th()*sample.sHe()*sample.sHe()/Math.sqrt(
                     (sample.He()*sample.He()*sample.sU()*sample.sU() + 
                     sample.U()*sample.U()*sample.sHe()*sample.sHe())*
                     (sample.He()*sample.He()*sample.sTh()*sample.sTh() + 
                     sample.Th()*sample.Th()*sample.sHe()*sample.sHe()));
        double[] xyz = {sample.U(),sample.Th(),sample.He()};
        double[] VWbar = Calculate.UThHe2vw(xyz);
        double sV = Math.sqrt((sample.sHe()/sample.He())*(sample.sHe()/sample.He()) + 
                (sample.sU()/sample.U())*(sample.sU()/sample.U()))/Math.log(10);
        double sW = Math.sqrt((sample.sHe()/sample.He())*(sample.sHe()/sample.He()) + 
                (sample.sTh()/sample.Th())*(sample.sTh()/sample.Th()))/Math.log(10);
        buildEllipse(n,numsigma,VWbar[0],sV,VWbar[1],sW,rho);
    }

    public Ellipse(Sample sample) throws Exception {
        this(10, 1, sample);
    }
    
    // this version of the constructor is for plotting average logratio compositions
    public Ellipse(int n, double numsigma, double V, double sV, double W, double sW, double covVW) throws Exception {
        double rho = covVW/(sV*sW);
        buildEllipse(n, numsigma, V, sV, W, sW, rho);
    }
    
    private void buildEllipse(int n, double m, double V, double sV, double W, double sW, double rho){
        this.n = n;        
        ellipse = new double[this.n][3]; // n rows, 4 columns: U, Th, He        
        double alpha = 0.5 * Math.atan(2 * rho * sV * sW / (sW*sW - sV*sV));
        double a = Math.sqrt(sV*sV * sW*sW * (1 - rho*rho) / (sV*sV * Math.cos(alpha)*Math.cos(alpha) - 
             2 * rho * sV * sW * Math.sin(alpha) * Math.cos(alpha) + 
             sW*sW * Math.sin(alpha)*Math.sin(alpha)));
        double b = Math.sqrt(sV*sV * sW*sW * (1 - rho*rho) / (sV*sV * Math.sin(alpha)*Math.sin(alpha) + 
                2 * rho * sV * sW * Math.sin(alpha) * Math.cos(alpha) + 
                sW*sW * Math.cos(alpha)*Math.cos(alpha)));
        double t;
        double[] VW = new double[2], xyz;
        for (int i=0; i<n; i++){
            t = i*2*Math.PI/(n-1);
            VW[0] = m * a * Math.cos(t) * Math.sin(alpha) + m * b * Math.sin(t) * Math.cos(alpha) + V;
            VW[1] = m * a * Math.cos(t) * Math.cos(alpha) - m * b * Math.sin(t) * Math.sin(alpha) + W;
            xyz = Calculate.vw2UThHe(VW);
            ellipse[i][0] = xyz[0];
            ellipse[i][1] = xyz[1];
            ellipse[i][2] = xyz[2];
        }        
    }
    
    public double[] getRow(int i) throws Exception {
        double[] p = {ellipse[i][0],ellipse[i][1],ellipse[i][2]};
        return p;
    }
    
    public int size(){
        return this.n;
    }

    private double[][] ellipse;
    private int n;
}