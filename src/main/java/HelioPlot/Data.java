package HelioPlot;

import Jama.Matrix;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.TreeMap;

public class Data implements Iterator, Iterable {

    public Data(){
        samplename = "";
        composition = new ArrayList[n];
        // initialise the data columns
        for (int i=0; i<n; i++){
            composition[i] = new ArrayList<Double>();
        }
    }

    /* boolean hasSm()
     * returns true if all Sm values are >0 and !=0
     */
    public boolean hasSm() throws Exception {
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            if (sample.Sm()<=0.0){
                return false;
            }
        }
        return true;
    }
    
    public boolean plotSm() throws Exception {
        return (!this.hasColour() & this.hasSm());
    }

    /* boolean hasColour()
     * returns true if there is any variation in Colour values
     */
    public boolean hasColour() throws Exception {
        Iterator i = this.iterator();
        Sample sample = (Sample)i.next();
        double C1 = sample.C();
        for (; i.hasNext();) {
            sample = (Sample)i.next();
            if (sample.C()!=C1){
                return true;
            }
        }
        return false;
    }    
    
    public void load(String filename) {
        String aLine;
        try {
            if (!filename.equals("")){
                System.out.println("Loading: " + filename);
                FileInputStream fin = new FileInputStream(filename);
                BufferedReader br = new BufferedReader(new InputStreamReader(fin));
                aLine = br.readLine();
                StringTokenizer stokenizer = new StringTokenizer(aLine, ",");
                samplename = stokenizer.nextToken();
                for (int i=0; i<n; i++){
                    composition[i].clear();
                }
                int numtokens;
                // extract columns of composition
                for (int i=0; (aLine = br.readLine()) != null ; i++) {
                    StringTokenizer st = new StringTokenizer(aLine, ",");
                    numtokens = (st.countTokens()<n) ? st.countTokens() : n;
                    for (int j=0; j<numtokens; j++){
                        try {
                            composition[j].add(Double.parseDouble(st.nextToken()));
                        } catch (NumberFormatException e){
                            composition[j].add(Data.NAN);
                        }
                    }
                    for(int j=numtokens; j<n; j++){
                        composition[j].add(Data.NAN);
                    }
                }
                br.close();
            }
        } catch (IOException ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }

    void write(String filepath) {
        String nl = System.getProperties().getProperty("line.separator");
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(filepath));
            out.write(samplename + nl);
            for (int i=0; i<composition[0].size(); i++){
                for (int j=0; j<n-1; j++){
                    out.write(composition[j].get(i) + ",");
                }
                out.write(composition[n-1].get(i) + nl);
            }
            out.close();
        } catch (IOException ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }

    /* double[][] averageCovVWX()
     * calculate the covariance matrix of the average logratio composition
     */
    public Matrix averageCovVWX() throws Exception{
        Sample sample;
        boolean doSm = this.hasSm();
        int m = doSm ? 3 : 2;
        Matrix sum = new Matrix(m,m), omega;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            omega = sample.cov(doSm).inverse();
            sum = sum.plusEquals(omega);
        }
        return sum.inverse();
    }

    public double MSWD() throws Exception {
        double mswd;
        int N = 0;
        Matrix vw = new Matrix(2,1), vwBar = this.mle2(),
               omega, S = new Matrix(1,1);
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            N++;
            sample = (Sample)i.next();
            vw.set(0, 0, sample.V());
            vw.set(1, 0, sample.W());
            omega = sample.cov(false).inverse();
            S = S.plus(vw.minus(vwBar).transpose().times(omega.times(vw.minus(vwBar))));
        }
        mswd = S.det()/(2*N-2);
        if (mswd > 0){
            return mswd;
        } else {
            return Data.NAN;
        }
    }

    public double MSWD(double tbar, boolean arithmetic) throws Exception {
        ArrayList<double[]> tst = arithmetic ? this.getAgeErrs() : this.getLogAgeErrs();
        double mswd = 0, t, st;
        int m = tst.size();
        for (int i=0; i<m; i++){
            t = tst.get(i)[0];
            st = tst.get(i)[1];
            if (arithmetic){
                mswd += (t-tbar)*(t-tbar)/(st*st)/(m-1);
            } else {
                mswd += (t-Math.log(tbar))*(t-Math.log(tbar))/(st*st)/(m-1);
            }
        }
        return mswd;
    }

    /* returns the 2-dimensional maximum likelihood mean assuming zero overdispersion
     */
    public Matrix mle2() throws Exception {
        Object[] v = this.getV().toArray(),
                 w = this.getW().toArray(),
                 varv = this.getVarV().toArray(),
                 varw = this.getVarW().toArray(),
                 covvw = this.getCovVW().toArray();
        return Newton.solveMu(new Matrix(2,2),v.length,varv,varw,covvw,v,w);
    }

    /* returns the 2-dimensional maximum likelihood mean assuming zero overdispersion
     */
    public Matrix mle3() throws Exception {
        Object[] v = this.getV().toArray(),
                 w = this.getW().toArray(),
                 x = this.getX().toArray(),
                 varv = this.getVarV().toArray(),
                 varw = this.getVarW().toArray(),
                 varx = this.getVarX().toArray(),
                 covvw = this.getCovVW().toArray(),
                 covvx = this.getCovVX().toArray(),
                 covwx = this.getCovWX().toArray();
        return Newton.solveMu(new Matrix(3,3),v.length,varv,varw,varx,covvw,covvx,covwx,v,w,x);
    }

    public double[] centralAgeErr() throws Exception {
        Matrix[] MuCov = Newton.solveMuCov(this,true);
        return Calculate.ageErr(MuCov);
    }

    public ArrayList<double[]> getAgeErrs() throws Exception {
        ArrayList<double[]> tst = new ArrayList<double[]>();
        Sample sample;
        for (Iterator i = iterator(); i.hasNext(); ){
            sample = (Sample)i.next();
            tst.add(Calculate.ageErr(sample));
        }
        return tst;
    }

    public ArrayList<double[]> getLogAgeErrs() throws Exception {
        ArrayList<double[]> tst = new ArrayList<double[]>();
        Sample sample;
        for (Iterator i = iterator(); i.hasNext(); ){
            sample = (Sample)i.next();
            tst.add(Calculate.logAgeErr(sample));
        }
        return tst;
    }

    double[] averageAgeErr(boolean arithmetic) throws Exception {
        double[] MuVar = Newton.solveMLEtv(this,arithmetic);
        double[] ae = {MuVar[0],Math.sqrt(MuVar[1])};
        return ae;
    }

    // returns the number of ages
    public int numAges(){
        int num = 0;
        for (Iterator i = this.iterator(); i.hasNext(); i.next()) {
            num++;
        }
        return num;
    }

    /* double[] CI(double p)
     * Bayesian Credibility Interval
     * input: p = percentile
     * output: [min max]
     */
    public double[] CI(double p) throws Exception {
        boolean doSm = this.hasSm();
        Matrix[] MuCov = Newton.solveMuCov(this);
        Matrix mu = MuCov[0],
               cov = MuCov[1],
               vwx = new Matrix(3,1);
        int m = 3, nV, nW, nX;
        double f, t, V, W, X, sV, sW, sX, dV, dW, dX,
               minV, maxV, minW, maxW, minX = NAN-1, maxX = NAN;
        TreeMap<Double,Double> pdf = new TreeMap(), cdf;
        double[] UThSmHe, out = new double[2];
        // end of variable declarations
        sV = Math.sqrt(cov.get(0, 0));
        sW = Math.sqrt(cov.get(1, 1));
        minV = mu.get(0, 0) - m*sV;
        maxV = mu.get(0, 0) + m*sV;
        minW = mu.get(1, 0) - m*sW;
        maxW = mu.get(1, 0) + m*sW;
        if (doSm & cov.det()<=0) { // to handle non positive definite covariances
            mu = mu.getMatrix(0,1,0,0);
            cov = cov.getMatrix(0,1,0,1);
            doSm = false;
        }
        if (doSm){
            nV = 15; nW = 15; nX = 10;
            sX = Math.sqrt(cov.get(2, 2));
            minX = mu.get(2, 0) - m*sX;
            maxX = mu.get(2, 0) + m*sX;
        } else {
            nV = 50; nW = 50; nX = 1;
        }
        dV = (maxV-minV)/nV;
        dW = (maxW-minW)/nW;
        dX = (maxX-minX)/nX;
        X = minX;
        while (X<maxX){
            X += dX;
            W = minW;
            while (W<maxW){
                W += dW;
                V = minV;
                while (V<maxV){
                    V += dV;
                    vwx.set(0, 0, V); vwx.set(1, 0, W); vwx.set(2, 0, X);
                    UThSmHe = Calculate.vwx2UThSmHe(vwx);
                    t = Calculate.age(UThSmHe[0], UThSmHe[1], UThSmHe[2], UThSmHe[3]);
                    if (doSm){
                        f = Stat.normalPDF(vwx, mu, cov);
                    } else {
                        f = Stat.normalPDF(vwx.getMatrix(0, 1, 0, 0), mu, cov);
                    }
                    if (pdf.containsKey(t)){
                        f += pdf.get(t);
                    }
                    pdf.put(t, f);
                }
            }
        }
        cdf = Stat.pdf2cdf(pdf);
        //Stat.printTree(cdf);
        out[0] = Stat.percentile(cdf,p/2);
        out[1] = Stat.percentile(cdf,1-p/2);
        return out;
    }

    void setSampleName(String samplename) {
        this.samplename = samplename;
    }

    public double getSmBar() throws Exception {
        double Sm = 0d;
        boolean doSm = this.hasSm();
        double[] UThSmHe;
        if (doSm){
            Matrix vwxBar = this.mle3();
            UThSmHe = Calculate.vwx2UThSmHe(vwxBar);
            Sm = UThSmHe[2];
        }
        return Sm;
    }

    public int size() {
        return composition[0].size();
    }

    public int numCols(){
        return n;
    }

    public void removeRow(int r) throws Exception {
        for (int i=0; i<n; i++){
            composition[i].remove(r);
        }
    }

    public void insertRow(int r) throws Exception {
        for (int i=0; i<n; i++){
            composition[i].add(r, NAN);
        }
    }

    public double get(int r, int c) throws Exception {
        return composition[c].get(r);
    }

    public void set(double val, int c, int r) throws Exception {
        // you won't go through this for loop if composition[c].size > r
        for (int i = composition[c].size(); i <= r; i++){
            for (int j=0; j<n; j++){
                composition[j].add(NAN);
            }
        }
        composition[c].set(r, val);
    }

    public double[] getRow(int r) {
        try {
            double[] row = new double[n];
            for (int i = 0 ; i<n; i++){
                row[i] = composition[i].get(r);
            }
            return row;
        } catch (Exception ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
            return new double[n];
        }
    }

    private double getU(int i){
        return composition[0].get(i);
    }

    private double getTh(int i){
        return composition[2].get(i);
    }

    private double getHe(int i){
        return composition[6].get(i);
    }

    public double getC(int i){
        return composition[8].get(i);
    }

    public ArrayList<Double> getSm() throws Exception {
        ArrayList<Double> Sm = new ArrayList<Double>();
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            Sm.add(sample.Sm());
        }
        return Sm;
    }

    public ArrayList<Double> getC() throws Exception {
        ArrayList<Double> C = new ArrayList<Double>();
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            C.add(sample.C());
        }
        return C;
    }
    
    

    public ArrayList<Double> getV() throws Exception {
        ArrayList<Double> V = new ArrayList<Double>();
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            V.add(sample.V());
        }
        return V;
    }

    public ArrayList<Double> getVarV() throws Exception {
        ArrayList<Double> vV = new ArrayList<Double>();
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            vV.add(sample.vV());
        }
        return vV;
    }

    public ArrayList<Double> getW() throws Exception {
        ArrayList<Double> W = new ArrayList<Double>();
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            W.add(sample.W());
        }
        return W;
    }

    public ArrayList<Double> getVarW() throws Exception {
        ArrayList<Double> vW = new ArrayList<Double>();
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            vW.add(sample.vV());
        }
        return vW;
    }

    public ArrayList<Double> getX() throws Exception {
        ArrayList<Double> X = new ArrayList<Double>();
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            X.add(sample.X());
        }
        return X;
    }

    public ArrayList<Double> getVarX() throws Exception {
        ArrayList<Double> vX = new ArrayList<Double>();
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            vX.add(sample.vX());
        }
        return vX;
    }

    public ArrayList<Double> getCovVW() throws Exception {
        ArrayList<Double> covVW = new ArrayList<Double>();
        Sample sample;
        for (Iterator i = this.iterator(); i.hasNext(); ) {
            sample = (Sample)i.next();
            covVW.add(sample.covVW());
        }
        return covVW;
    }

    public ArrayList<Double> getCovVX() throws Exception{
        return this.getCovVW();
    }

    public ArrayList<Double> getCovWX() throws Exception{
        return this.getCovVW();
    }

    @Override
    public boolean hasNext() {
        return ii < this.size();
    }

    @Override
    public Object next() {
        Sample sample = new Sample(this.getRow(ii));
        // put iterator in right position for next request
        do {
            ii++;
        } while (hasNext() && (this.getU(ii) == NAN || this.getTh(ii) == NAN || this.getHe(ii) == NAN));
        return sample;
    }

    int getRowNumber(int samplenumber) throws Exception {
        int sn = 0;
        for (Iterator i = iterator(); i.hasNext(); i.next()) {
            if (sn >= samplenumber){
                break;
            } else {
                sn++;
            }
        }
        return ii;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Data:Iterator:remove() not supported");
    }

    @Override
    public Iterator iterator() {
        ii = 0;
        // if the data are empty, return immediately
        if (this.size()==0){return this;}
        // if the first row is empty, find the first non-empty row
        if (this.getU(ii)<=0 || this.getTh(ii)<=0 || this.getHe(ii)<=0){
            do {
            ii++;
            } while (hasNext() && (this.getU(ii)<=0 || this.getTh(ii)<=0 || this.getHe(ii)<=0));
            return this;
        } else {
            return this;
        }
    }

    protected String samplename;
    protected ArrayList<Double>[] composition;
    static double NAN = -999.9;
    static final int n = 9;
    protected int ii;

}