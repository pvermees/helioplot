package HelioPlot;

import Jama.Matrix;
import java.util.ArrayList;
import java.util.Iterator;

public class Newton {

static Matrix[] solveMuCov(Data data) throws Exception {
    return solveMuCov(data, false);
}

/* iteratively solve for maximum likelihood compositional mean and covariance
 */
static Matrix[] solveMuCov(Data data, boolean printmessage) throws Exception {
    Matrix[] MuCov = new Matrix[2];
    String msg;
    Object[] v = data.getV().toArray(),
             w = data.getW().toArray(),
             x = data.getX().toArray(),
             varv = data.getVarV().toArray(),
             varw = data.getVarW().toArray(),
             varx = data.getVarX().toArray(),
             covvw = data.getCovVW().toArray(),
             covvx = data.getCovVX().toArray(),
             covwx = data.getCovWX().toArray();
    boolean doSm = data.hasSm();
    int n = data.getV().size();
    try {
        Matrix[] MuXi = solveMuXi(v,w,x,varv,varw,varx,covvw,covvx,covwx,doSm,n);
        MuCov[0] = MuXi[0];
        MuCov[1] = getCov(MuXi[1], data);
        msg = "Central age: internal + external error.";
        throw new Exception();
    } catch (Exception e){
        double mswd = data.MSWD();
        if (mswd < 3){
            // assume zero overdispersion
            int m = doSm ? 3 : 2;
            Matrix Xi = new Matrix(m,m); // initialise to zero
            MuCov[0] = solveMu(Xi,v,w,x,varv,varw,varx,covvw,covvx,covwx,doSm,n);
            MuCov[1] = getCov(Xi,data);
            msg = "Central age: internal error only.";
        } else {
            // assume infinite precision
            Matrix Xi = initXi(v,w,x,varv,varw,varx,covvw,covvx,covwx,doSm,n);
            MuCov[0] = initMu(Xi,v,w,x,doSm,n);
            MuCov[1] = Xi.inverse().times(data.numAges()).inverse();
            msg = "Central age: external error only.";
        }
    }
    if (printmessage){System.out.println(msg);}
    return MuCov;
}

/* iteratively solve for maximum likelihood age mean and variance
 */
static double[] solveMLEtv(Data data, boolean arithmetic) throws Exception {
    ArrayList<double[]> tst = arithmetic ? data.getAgeErrs() : data.getLogAgeErrs();
    double[] MuXi = solveMuXi(tst),
             MuVar = new double[2];
    double Var = gettv(tst,MuXi[1]);
    if (arithmetic){
        MuVar[0] = MuXi[0];
        MuVar[1] = Var;
    } else {
        MuVar[0] = Math.exp(MuXi[0]);
        MuVar[1] = Var*MuVar[0]*MuVar[0];
    }
    return MuVar;
}

/* calculates age average and overdispersion
 */
static double[] solveMuXi(ArrayList<double[]> tst) throws Exception {
    double[] MuXi = new double[2];
    double Mu = 0d, Xi = 0d, dMu = 0d;
    for (int i=0; i<n2; i++){
        MuXi[0] = Mu;
        MuXi[1] = Xi;
        Mu = solveMu(tst,Xi);
        Xi = solveXi(tst,Mu);
        dMu = Mu - MuXi[0];
        if (doBreak(dMu,MuXi[0])){
            break;
        }
    }
    return MuXi;
}

static double solveMu(ArrayList<double[]> tst, double Xi) throws Exception {
    double ti, si, num = 0d, denom = 0d;
    for (int i=0; i<tst.size(); i++){
        ti = tst.get(i)[0];
        si = tst.get(i)[1];
        num += ti/(si*si + Xi);
        denom += 1/(si*si + Xi);
    }
    return num/denom;
}

static double solveXi(ArrayList<double[]> tst, double mu) throws Exception {
    double xi = 0d, dxi;
    double[] fdf = fdf = getfdf(tst,mu,xi);
    // exit with zero if the first iteration sends you to negative xi
    if (fdf[0]/fdf[1] > 0){return xi;}
    for (int i=0; i<n2; i++){
        fdf = getfdf(tst,mu,xi);
        dxi = fdf[0]/fdf[1];
        xi -= dxi;
        // xi must be positive
        if (xi<0 || Math.sqrt(xi)/mu > 100){
            return 0;
        }
        if (doBreak(dxi,xi)){
            break;
        }
    }
    return xi;
}

static double[] getfdf(ArrayList<double[]> tst, double mu, double xi){
    double f = 0d, df = 0d, a, a2, b, b2, b3, t, st;
    for (int i=0; i<tst.size(); i++){
        t = tst.get(i)[0];
        st = tst.get(i)[1];
        a = t-mu;
        a2 = a*a;
        b = st*st+xi;
        b2 = b*b;
        b3 = b2*b;
        f += (a2/b2 - 1/b);
        df += (1/b2 - 2*a2/b3);
    }
    double[] fdf = {f,df};
    return fdf;
}

/* Calculates compositional average and overdispersion
 */
static Matrix[] solveMuXi(Object[] v, Object[] w, Object[] x,
                          Object[] varv, Object[] varw, Object[] varx,
                          Object[] covvw, Object[] covvx, Object[] covwx,
                          boolean doSm, int n) throws Exception {
    int m = doSm ? 3 : 2;
    Matrix[] MuXi = {new Matrix(m,1), new Matrix(m,m)};
    Matrix Mu = new Matrix(m,1), Xi = new Matrix(m,m);
    Matrix dXi;
    for (int i=0; ; i++){
        MuXi[0] = Mu;
        MuXi[1] = Xi;
        if (doSm){
            Xi = solveXi(MuXi,n,varv,varw,varx,covvw,covvx,covwx,v,w,x);
            Mu = solveMu(MuXi[1],n,varv,varw,varx,covvw,covvx,covwx,v,w,x);
        } else {
            Xi = solveXi(MuXi,n,varv,varw,covvw,v,w);
            Mu = solveMu(MuXi[1],n,varv,varw,covvw,v,w);
        }
        dXi = Xi.minus(MuXi[1]);
        if (doBreak(dXi,MuXi[1])){break;}
        if (i==n1){throw new Exception();}
    }
    return MuXi;
}

private static double gettv(ArrayList<double[]> tst, double Xi) throws Exception {
    double denom = 0d, st;
    for (int i=0; i<tst.size(); i++){
        st = tst.get(i)[1];
        denom += 1/(st*st + Xi);
    }
    return 1/denom;
}

static Matrix getCov(Matrix Xi, Data data) throws Exception {
    boolean doSm = data.hasSm();
    int m = doSm ? 3 : 2;
    Matrix EplusXi, cov, foo = new Matrix(m,m);
    for (Iterator<Sample> i = data.iterator(); i.hasNext(); ) {
        Sample sample = i.next();
        cov = sample.cov(doSm);
        EplusXi = cov.plus(Xi);
        foo.plusEquals(EplusXi.inverse());
    }
    return foo.inverse();
}

static Matrix solveMu(Matrix Xi, Object[] v, Object[] w, Object[] x,
                      Object[] varv, Object[] varw, Object[] varx,
                      Object[] covvw, Object[] covvx, Object[] covwx,
                      boolean doSm, int n) throws Exception{
    if (doSm){
        return solveMu(Xi,n,varv,varw,varx,covvw,covvx,covwx,v,w,x);
    } else {
        return solveMu(Xi,n,varv,varw,covvw,v,w);
    }
}

static Double[] zeros(int n){
    Double[] zeros = new Double[n];
    for (int i=0;i<n;i++){zeros[i] = Double.valueOf(0d);}
    return zeros;
}

static Matrix solveMu(Matrix Xi, int n, Object[] v, Object[] w) throws Exception {
    Double[] zeros = zeros(n);
    return solveMu(Xi,n,zeros,zeros,zeros,v,w);
}

/* find Mu given Xi and the analytical precisions */
static Matrix solveMu(Matrix Xi, int n, Object[] varv, Object[] varw,
                      Object[] covvw, Object[] v, Object[] w) throws Exception {
    Matrix A = new Matrix(2,2), E = new Matrix(2,2);
    Matrix foo, bar, B = new Matrix(2,1), X = new Matrix(2,1);
    for (int i=0; i<n; i++){
        E.set(0,0,((Double)varv[i]).doubleValue());
        E.set(0,1,((Double)covvw[i]).doubleValue());
        E.set(1,0,((Double)covvw[i]).doubleValue());
        E.set(1,1,((Double)varw[i]).doubleValue());
        X.set(0,0,((Double)v[i]).doubleValue());
        X.set(1,0,((Double)w[i]).doubleValue());
        foo = E.plus(Xi).inverse();
        bar = foo.times(X);
        A.plusEquals(foo);
        B.plusEquals(bar);
    }
    return A.solve(B);
}

static Matrix solveXi(Matrix[] MuXi, int n, Object[] varv, Object[] varw,
                     Object[] covvw, Object[] v, Object[] w) throws Exception {
    Matrix Mu = MuXi[0], Xi = new Matrix(2,2);
    Matrix F = new Matrix(3,1), J = new Matrix(3,3),
            S = new Matrix(3,1), dS = new Matrix(3,1);
    for (int i=0; i<n2; i++){
        F = F3x3(S.get(0,0),S.get(2,0),S.get(1,0),Mu.get(0,0),Mu.get(1,0),n,v,w,varv,varw,covvw);
        J = J3x3(S.get(0,0),S.get(2,0),S.get(1,0),Mu.get(0,0),Mu.get(1,0),n,v,w,varv,varw,covvw);
        dS = J.solve(F).times(-1);
        S.plusEquals(dS);
        if (doBreak(dS,S)){
            break;
        }
    }
    Xi.set(0,0,S.get(0,0)); Xi.set(1,1,S.get(1,0));
    Xi.set(0,1,S.get(2,0)); Xi.set(1,0,S.get(2,0));
    return Xi;
}

static Matrix solveMu(Matrix Xi, int n, Object[] v, Object[] w, Object[] x) throws Exception {
    Double[] zeros = zeros(n);
    return solveMu(Xi,n,zeros,zeros,zeros,zeros,zeros,zeros,v,w,x);
}

static Matrix solveMu(Matrix Xi, int n, Object[] varv, Object[] varw, Object[] varx,
                      Object[] covvw, Object[] covvx, Object[] covwx,
                      Object[] v, Object[] w, Object[] x) throws Exception {
    Matrix A = new Matrix(3,3), E = new Matrix(3,3);
    Matrix foo, bar, B = new Matrix(3,1), X = new Matrix(3,1);
    for (int i=0; i<n; i++){
        E.set(0,0,((Double)varv[i]).doubleValue());
        E.set(0,1,((Double)covvw[i]).doubleValue());
        E.set(0,2,((Double)covvx[i]).doubleValue());
        E.set(1,0,((Double)covvw[i]).doubleValue());
        E.set(1,1,((Double)varw[i]).doubleValue());
        E.set(1,2,((Double)covwx[i]).doubleValue());
        E.set(2,0,((Double)covvx[i]).doubleValue());
        E.set(2,1,((Double)covwx[i]).doubleValue());
        E.set(2,2,((Double)varx[i]).doubleValue());
        X.set(0,0,((Double)v[i]).doubleValue());
        X.set(1,0,((Double)w[i]).doubleValue());
        X.set(2,0,((Double)x[i]).doubleValue());
        foo = E.plus(Xi).inverse();
        bar = foo.times(X);
        A.plusEquals(foo);
        B.plusEquals(bar);
    }
    return A.solve(B);
}

static Matrix solveXi(Matrix[] MuXi, int n,
        Object[] varv, Object[] varw, Object[] varx,
        Object[] covvw, Object[] covvx, Object[] covwx,
        Object[] v, Object[] w, Object[] x) throws Exception {
    Matrix Mu = MuXi[0], Xi = new Matrix(3,3);
    Matrix F = new Matrix(6,1), J = new Matrix(6,6),
           S = new Matrix(6,1), dS = new Matrix(6,1);
    for (int i=0; i<n2; i++){
        F = F6x6(S.get(0,0),S.get(1,0),S.get(2,0),S.get(3,0),S.get(4,0),S.get(5,0),
                 Mu.get(0,0),Mu.get(1,0),Mu.get(2,0),n,v,w,x,varv,varw,varx,covvw,covvx,covwx);
        J = J6x6(S.get(0,0),S.get(1,0),S.get(2,0),S.get(3,0),S.get(4,0),S.get(5,0),
                 Mu.get(0,0),Mu.get(1,0),Mu.get(2,0),n,v,w,x,varv,varw,varx,covvw,covvx,covwx);
        dS = J.solve(F).times(-1);
        S.plusEquals(dS);
        if (doBreak(dS,S)){
            break;
        }
    }
    Xi.set(0,0,S.get(0,0)); Xi.set(0,1,S.get(3,0)); Xi.set(0,2,S.get(4,0));
    Xi.set(1,0,S.get(3,0)); Xi.set(1,1,S.get(1,0)); Xi.set(1,2,S.get(5,0));
    Xi.set(2,0,S.get(4,0)); Xi.set(2,1,S.get(5,0)); Xi.set(2,2,S.get(2,0));
    return Xi;
}

static private boolean doBreak(Matrix dS, Matrix S) throws Exception {
    return doBreak(dS.norm1(),S.norm1());
}

static private boolean doBreak(double ds, double s) throws Exception {
    boolean out = (s!=0 && Math.abs(ds/s)<1e-5);
    return out;
}

// initialize Mu ignoring the analytical precision
static Matrix initMu(Matrix Xi, Object[] v, Object[] w, Object[] x,
                     boolean doSm, int n) throws Exception {
    int m = doSm ? 3 : 2;
    if (m==2){
        return solveMu(Xi,n,v,w);
    } else {
        return solveMu(Xi,n,v,w,x);
    }
}

/* initialize Xi ignoring the analytical error
   (if more than data contains more than one measurement
   otherwise use measurement covariance) */
static Matrix initXi(Object[] v, Object[] w, Object[] x,
                    Object[] varv, Object[] varw, Object[] varx,
                    Object[] covvw, Object[] covvx, Object[] covwx,
                    boolean doSm, int n) throws Exception {
    if (doSm){
        return initXi(n,varv,varw,varx,covvw,covvx,covwx,v,w,x);
    } else {
        return initXi(n,varv,varw,covvw,v,w);
    }
}

// set Xi equal to the covariance matrix of the logratio data
static Matrix initXi(int n, Object[] varv, Object[] varw,
                   Object[] covvw, Object[] v, Object[] w) throws Exception {
    Matrix Xi = new Matrix(2,2);
    if (n==1){
        Xi.set(0,0,(Double)varv[0]);
        Xi.set(1,1,(Double)varw[0]);
        Xi.set(0,1,(Double)covvw[0]);
        Xi.set(1,0,Xi.get(0,1));
    } else {
        Xi.set(0,0,Stat.var(n,v));
        Xi.set(1,1,Stat.var(n,w));
        Xi.set(0,1,Stat.cov(n,v,w));
        Xi.set(1,0,Xi.get(0,1));
    }
    return Xi;
}

// set Xi equal to the covariance matrix of the logratio data
static Matrix initXi(int n, Object[] varv, Object[] varw, Object[] varx,
        Object[] covvw, Object[] covvx, Object[] covwx,
        Object[] v, Object[] w, Object[] x) throws Exception {
    Matrix Xi = new Matrix(3,3);
    if (n==1){
        Xi.set(0,0,(Double)varv[0]);
        Xi.set(1,1,(Double)varw[0]);
        Xi.set(2,2,(Double)varx[0]);
        Xi.set(0,1,(Double)covvw[0]);
        Xi.set(0,2,(Double)covvx[0]);
        Xi.set(1,2,(Double)covwx[0]);
        Xi.set(1,0,Xi.get(0,1));
        Xi.set(2,0,Xi.get(0,2));
        Xi.set(2,1,Xi.get(2,1));
    } else {
        Xi.set(0,0,Stat.var(n,v));
        Xi.set(1,1,Stat.var(n,w));
        Xi.set(2,2,Stat.var(n,x));
        Xi.set(0,1,Stat.cov(n,v,w));
        Xi.set(0,2,Stat.cov(n,v,x));
        Xi.set(1,2,Stat.cov(n,w,x));
        Xi.set(1,0,Xi.get(0,1));
        Xi.set(2,0,Xi.get(0,2));
        Xi.set(2,1,Xi.get(2,1));
    }
    return Xi;
}

static Matrix F3x3(double S11, double S22, double S12, double muv, double muw, int n,
        Object[] v, Object[] w, Object[] varv, Object[] varw, Object[] covvw) throws Exception {
    Matrix F = new Matrix(3,1);
    for (int i=0; i<n; i++){
        F.plusEquals(getF3x3(S11,S22,S12,muv,muw,
                ((Double)v[i]).doubleValue(),
                ((Double)w[i]).doubleValue(),
                ((Double)varv[i]).doubleValue(),
                ((Double)varw[i]).doubleValue(),
                ((Double)covvw[i]).doubleValue()));
    }
    return F;
}

static Matrix getF3x3(double S11, double S22, double S12, double muv, double muw,
        double v, double w, double varv, double varw, double covvw) throws Exception {
    Matrix F = new Matrix(3,1);
    double A = (varv*varw+varv*S22+S11*varw+S11*S22-(covvw*covvw)-2*covvw*S12-(S12*S12));
    F.set(0,0,((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)*(varw+S22)/A-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)*(covvw+S12)/A-(varw+S22)/A);
    F.set(1,0,-(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(v-muv)*(covvw+S12)/A+(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(w-muw)*(varv+S11)/A-(varv+S11)/A);
    F.set(2,0,-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)*(covvw+S12)/A+((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)*(varv+S11)/A+(covvw+S12)/A);
    return F;
}

static Matrix F6x6(double S11, double S22, double S33, double S12, double S13, double S23,
        double muv, double muw, double mux, int n,
        Object[] v, Object[] w, Object[] x, Object[] varv, Object[] varw, Object[] varx,
        Object[] covvw, Object[] covvx, Object[] covwx) throws Exception {
    Matrix F = new Matrix(6,1);
    for (int i=0; i<n; i++){
        F.plusEquals(getF6x6(S11,S22,S33,S12,S13,S23,muv,muw,mux,
                             ((Double)v[i]).doubleValue(),((Double)w[i]).doubleValue(),
                             ((Double)x[i]).doubleValue(),((Double)varv[i]).doubleValue(),
                             ((Double)varw[i]).doubleValue(),((Double)varx[i]).doubleValue(),
                             ((Double)covvw[i]).doubleValue(),((Double)covvx[i]).doubleValue(),
                             ((Double)covwx[i]).doubleValue()));
    }
    return F;
}

static Matrix getF6x6(double S11, double S22, double S33, double S12, double S13, double S23,
                double muv, double muw, double mux, double v, double w, double x,
                double varv, double varw, double varx, double covvw, double covvx, double covwx)
                throws Exception {
    Matrix f = new Matrix(6,1);
    double A,B,C,D,E,F,G;
    A = (varw*varx+varw*S33+S22*varx+S22*S33-(covwx*covwx)-2*covwx*S23-(S23*S23));
    B = (S11*varw*S33-2*covvw*S12*varx+2*S12*S13*S23+S11*varw*varx-2*covvx*S13*varw-(covvw*covvw)*varx-(covvw*covvw)*S33+varv*varw*varx+varv*varw*S33+varv*S22*varx+varv*S22*S33-2*varv*covwx*S23+S11*S22*varx+S11*S22*S33-2*S11*covwx*S23-2*covvw*S12*S33+2*covvw*covvx*covwx+2*covvw*covvx*S23+2*covvw*S13*covwx+2*covvw*S13*S23+2*S12*covvx*covwx+2*S12*covvx*S23+2*S12*S13*covwx-2*covvx*S13*S22-varv*(covwx*covwx)-varv*(S23*S23)-S11*(covwx*covwx)-S11*(S23*S23)-(S12*S12)*varx-(S12*S12)*S33-(covvx*covvx)*varw-(covvx*covvx)*S22-(S13*S13)*varw-(S13*S13)*S22);
    C = (covvw*varx+covvw*S33+S12*varx+S12*S33-covvx*covwx-covvx*S23-S13*covwx-S13*S23);
    D = (covvw*covwx+covvw*S23+S12*covwx+S12*S23-covvx*varw-covvx*S22-S13*varw-S13*S22);
    E = (varv*varx+varv*S33+S11*varx+S11*S33-(covvx*covvx)-2*covvx*S13-(S13*S13));
    F = (varv*covwx+varv*S23+S11*covwx+S11*S23-covvw*covvx-covvw*S13-S12*covvx-S12*S13);
    G = (varv*varw+varv*S22+S11*varw+S11*S22-(covvw*covvw)-2*covvw*S12-(S12*S12));
    f.set(0,0,(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*A/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*C/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*D/B-A/B);
    f.set(1,0,-(-C/B*(v-muv)+E/B*(w-muw)-F/B*(x-mux))*(v-muv)*C/B+(-C/B*(v-muv)+E/B*(w-muw)-F/B*(x-mux))*(w-muw)*E/B-(-C/B*(v-muv)+E/B*(w-muw)-F/B*(x-mux))*(x-mux)*F/B-E/B);
    f.set(2,0,(D/B*(v-muv)-F/B*(w-muw)+G/B*(x-mux))*(v-muv)*D/B-(D/B*(v-muv)-F/B*(w-muw)+G/B*(x-mux))*(w-muw)*F/B+(D/B*(v-muv)-F/B*(w-muw)+G/B*(x-mux))*(x-mux)*G/B-G/B);
    f.set(3,0,-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*C/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*E/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*F/B+C/B);
    f.set(4,0,(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*D/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*F/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*G/B-D/B);
    f.set(5,0,(-C/B*(v-muv)+E/B*(w-muw)-F/B*(x-mux))*(v-muv)*D/B-(-C/B*(v-muv)+E/B*(w-muw)-F/B*(x-mux))*(w-muw)*F/B+(-C/B*(v-muv)+E/B*(w-muw)-F/B*(x-mux))*(x-mux)*G/B+F/B);
    return f;
}

static Matrix J3x3(double S11, double S22, double S12, double muv, double muw, int n,
        Object[] v, Object[] w, Object[] varv, Object[] varw, Object[] covvw) throws Exception {
    Matrix J = new Matrix(3,3);
    for (int i=0; i<n; i++){
        J.plusEquals(getJ3x3(S11,S22,S12,muv,muw,
                ((Double)v[i]).doubleValue(),
                ((Double)w[i]).doubleValue(),
                ((Double)varv[i]).doubleValue(),
                ((Double)varw[i]).doubleValue(),
                ((Double)covvw[i]).doubleValue()));
    }
    return J;
}

static Matrix J6x6(double S11, double S22, double S33, double S12, double S13, double S23,
        double muv, double muw, double mux, int n,
        Object[] v, Object[] w, Object[] x, Object[] varv, Object[] varw, Object[] varx,
        Object[] covvw, Object[] covvx, Object[] covwx) throws Exception {
    Matrix J = new Matrix(6,6);
    for (int i=0; i<n; i++){
        J.plusEquals(getJ6x6(S11,S22,S33,S12,S13,S23,muv,muw,mux,
                ((Double)v[i]).doubleValue(),((Double)w[i]).doubleValue(),
                             ((Double)x[i]).doubleValue(),((Double)varv[i]).doubleValue(),
                             ((Double)varw[i]).doubleValue(),((Double)varx[i]).doubleValue(),
                             ((Double)covvw[i]).doubleValue(),((Double)covvx[i]).doubleValue(),
                             ((Double)covwx[i]).doubleValue()));
    }
    return J;
}

static Matrix getJ3x3(double S11, double S22, double S12, double muv, double muw,
        double v, double w, double varv, double varw, double covvw) throws Exception {
    Matrix J = new Matrix(3,3);
    double A = (varv*varw+varv*S22+S11*varw+S11*S22-covvw*covvw-2*covvw*S12-S12*S12);
    J.set(0,0,(-((varw+S22)*(varw+S22))/(A*A)*(v-muv)+(covvw+S12)/(A*A)*(w-muw)*(varw+S22))*(v-muv)*(varw+S22)/A-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)*((varw+S22)*(varw+S22))/(A*A)-(-((varw+S22)*(varw+S22))/(A*A)*(v-muv)+(covvw+S12)/(A*A)*(w-muw)*(varw+S22))*(w-muw)*(covvw+S12)/A+((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)*(covvw+S12)/(A*A)*(varw+S22)+((varw+S22)*(varw+S22))/(A*A));
    J.set(0,1,(1/A*(v-muv)-(varw+S22)/(A*A)*(v-muv)*(varv+S11)+(covvw+S12)/(A*A)*(w-muw)*(varv+S11))*(v-muv)*(varw+S22)/A+((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)/A-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)*(varw+S22)/(A*A)*(varv+S11)-(1/A*(v-muv)-(varw+S22)/(A*A)*(v-muv)*(varv+S11)+(covvw+S12)/(A*A)*(w-muw)*(varv+S11))*(w-muw)*(covvw+S12)/A+((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)*(covvw+S12)/(A*A)*(varv+S11)-1/A+(varw+S22)/(A*A)*(varv+S11));
    J.set(0,2,(-(varw+S22)/(A*A)*(v-muv)*(-2*covvw-2*S12)-1/A*(w-muw)+(covvw+S12)/(A*A)*(w-muw)*(-2*covvw-2*S12))*(v-muv)*(varw+S22)/A-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)*(varw+S22)/(A*A)*(-2*covvw-2*S12)-(-(varw+S22)/(A*A)*(v-muv)*(-2*covvw-2*S12)-1/A*(w-muw)+(covvw+S12)/(A*A)*(w-muw)*(-2*covvw-2*S12))*(w-muw)*(covvw+S12)/A-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)/A+((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)*(covvw+S12)/(A*A)*(-2*covvw-2*S12)+(varw+S22)/(A*A)*(-2*covvw-2*S12));
    J.set(1,0,-((covvw+S12)/(A*A)*(v-muv)*(varw+S22)+1/A*(w-muw)-(varv+S11)/(A*A)*(w-muw)*(varw+S22))*(v-muv)*(covvw+S12)/A+(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(v-muv)*(covvw+S12)/(A*A)*(varw+S22)+((covvw+S12)/(A*A)*(v-muv)*(varw+S22)+1/A*(w-muw)-(varv+S11)/(A*A)*(w-muw)*(varw+S22))*(w-muw)*(varv+S11)/A+(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(w-muw)/A-(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(w-muw)*(varv+S11)/(A*A)*(varw+S22)-1/A+(varw+S22)/(A*A)*(varv+S11));
    J.set(1,1,-((covvw+S12)/(A*A)*(v-muv)*(varv+S11)-((varv+S11)*(varv+S11))/(A*A)*(w-muw))*(v-muv)*(covvw+S12)/A+(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(v-muv)*(covvw+S12)/(A*A)*(varv+S11)+((covvw+S12)/(A*A)*(v-muv)*(varv+S11)-((varv+S11)*(varv+S11))/(A*A)*(w-muw))*(w-muw)*(varv+S11)/A-(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(w-muw)*((varv+S11)*(varv+S11))/(A*A)+((varv+S11)*(varv+S11))/(A*A));
    J.set(1,2,-(-1/A*(v-muv)+(covvw+S12)/(A*A)*(v-muv)*(-2*covvw-2*S12)-(varv+S11)/(A*A)*(w-muw)*(-2*covvw-2*S12))*(v-muv)*(covvw+S12)/A-(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(v-muv)/A+(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(v-muv)*(covvw+S12)/(A*A)*(-2*covvw-2*S12)+(-1/A*(v-muv)+(covvw+S12)/(A*A)*(v-muv)*(-2*covvw-2*S12)-(varv+S11)/(A*A)*(w-muw)*(-2*covvw-2*S12))*(w-muw)*(varv+S11)/A-(-(covvw+S12)/A*(v-muv)+(varv+S11)/A*(w-muw))*(w-muw)*(varv+S11)/(A*A)*(-2*covvw-2*S12)+(varv+S11)/(A*A)*(-2*covvw-2*S12));
    J.set(2,0,-(-((varw+S22)*(varw+S22))/(A*A)*(v-muv)+(covvw+S12)/(A*A)*(w-muw)*(varw+S22))*(v-muv)*(covvw+S12)/A+((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)*(covvw+S12)/(A*A)*(varw+S22)+(-((varw+S22)*(varw+S22))/(A*A)*(v-muv)+(covvw+S12)/(A*A)*(w-muw)*(varw+S22))*(w-muw)*(varv+S11)/A+((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)/A-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)*(varv+S11)/(A*A)*(varw+S22)-(covvw+S12)/(A*A)*(varw+S22));
    J.set(2,1,-(1/A*(v-muv)-(varw+S22)/(A*A)*(v-muv)*(varv+S11)+(covvw+S12)/(A*A)*(w-muw)*(varv+S11))*(v-muv)*(covvw+S12)/A+((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)*(covvw+S12)/(A*A)*(varv+S11)+(1/A*(v-muv)-(varw+S22)/(A*A)*(v-muv)*(varv+S11)+(covvw+S12)/(A*A)*(w-muw)*(varv+S11))*(w-muw)*(varv+S11)/A-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)*((varv+S11)*(varv+S11))/(A*A)-(covvw+S12)/(A*A)*(varv+S11));
    J.set(2,2,-(-(varw+S22)/(A*A)*(v-muv)*(-2*covvw-2*S12)-1/A*(w-muw)+(covvw+S12)/(A*A)*(w-muw)*(-2*covvw-2*S12))*(v-muv)*(covvw+S12)/A-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)/A+((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(v-muv)*(covvw+S12)/(A*A)*(-2*covvw-2*S12)+(-(varw+S22)/(A*A)*(v-muv)*(-2*covvw-2*S12)-1/A*(w-muw)+(covvw+S12)/(A*A)*(w-muw)*(-2*covvw-2*S12))*(w-muw)*(varv+S11)/A-((varw+S22)/A*(v-muv)-(covvw+S12)/A*(w-muw))*(w-muw)*(varv+S11)/(A*A)*(-2*covvw-2*S12)+1/A-(covvw+S12)/(A*A)*(-2*covvw-2*S12));
    return J;
}

static Matrix getJ6x6(double S11, double S22, double S33, double S12, double S13, double S23,
                double muv, double muw, double mux, double v, double w, double x,
                double varv, double varw, double varx, double covvw, double covvx, double covwx)
                throws Exception {
    Matrix j = new Matrix(6,6);
    double A,B,C,D,E,F,G,H,I,J;
    A = (varw*varx+varw*S33+S22*varx+S22*S33-(covwx*covwx)-2*covwx*S23-(S23*S23));
    B = (S11*varw*S33-2*covvw*S12*varx+2*S12*S13*S23+S11*varw*varx-2*covvx*S13*varw-(covvw*covvw)*varx-(covvw*covvw)*S33+varv*varw*varx+varv*varw*S33+varv*S22*varx+varv*S22*S33-2*varv*covwx*S23+S11*S22*varx+S11*S22*S33-2*S11*covwx*S23-2*covvw*S12*S33+2*covvw*covvx*covwx+2*covvw*covvx*S23+2*covvw*S13*covwx+2*covvw*S13*S23+2*S12*covvx*covwx+2*S12*covvx*S23+2*S12*S13*covwx-2*covvx*S13*S22-varv*(covwx*covwx)-varv*(S23*S23)-S11*(covwx*covwx)-S11*(S23*S23)-(S12*S12)*varx-(S12*S12)*S33-(covvx*covvx)*varw-(covvx*covvx)*S22-(S13*S13)*varw-(S13*S13)*S22);
    C = (covvw*varx+covvw*S33+S12*varx+S12*S33-covvx*covwx-covvx*S23-S13*covwx-S13*S23);
    D = (covvw*covwx+covvw*S23+S12*covwx+S12*S23-covvx*varw-covvx*S22-S13*varw-S13*S22);
    E = (varv*varx+varv*S33+S11*varx+S11*S33-(covvx*covvx)-2*covvx*S13-(S13*S13));
    F = (varv*varw+varv*S22+S11*varw+S11*S22-(covvw*covvw)-2*covvw*S12-(S12*S12));
    G = (-2*covvw*varx+2*S13*S23-2*covvw*S33+2*covvx*covwx+2*covvx*S23+2*S13*covwx-2*S12*varx-2*S12*S33);
    H = (2*S12*S23-2*covvx*varw+2*covvw*covwx+2*covvw*S23+2*S12*covwx-2*covvx*S22-2*S13*varw-2*S13*S22);
    I = (2*S12*S13-2*varv*covwx-2*S11*covwx+2*covvw*covvx+2*covvw*S13+2*S12*covvx-2*varv*S23-2*S11*S23);
    J = (varv*covwx+varv*S23+S11*covwx+S11*S23-covvw*covvx-covvw*S13-S12*covvx-S12*S13);
    j.set(0,0,(-(A*A)/(B*B)*(v-muv)+C/(B*B)*(w-muw)*A-D/(B*B)*(x-mux)*A)*(v-muv)*A/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(A*A)/(B*B)-(-(A*A)/(B*B)*(v-muv)+C/(B*B)*(w-muw)*A-D/(B*B)*(x-mux)*A)*(w-muw)*C/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*C/(B*B)*A+(-(A*A)/(B*B)*(v-muv)+C/(B*B)*(w-muw)*A-D/(B*B)*(x-mux)*A)*(x-mux)*D/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*D/(B*B)*A+(A*A)/(B*B));
    j.set(0,1,((varx+S33)/B*(v-muv)-A/(B*B)*(v-muv)*E+C/(B*B)*(w-muw)*E+(-covvx-S13)/B*(x-mux)-D/(B*B)*(x-mux)*E)*(v-muv)*A/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(varx+S33)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*A/(B*B)*E-((varx+S33)/B*(v-muv)-A/(B*B)*(v-muv)*E+C/(B*B)*(w-muw)*E+(-covvx-S13)/B*(x-mux)-D/(B*B)*(x-mux)*E)*(w-muw)*C/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*C/(B*B)*E+((varx+S33)/B*(v-muv)-A/(B*B)*(v-muv)*E+C/(B*B)*(w-muw)*E+(-covvx-S13)/B*(x-mux)-D/(B*B)*(x-mux)*E)*(x-mux)*D/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(-covvx-S13)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*D/(B*B)*E-(varx+S33)/B+A/(B*B)*E);
    j.set(0,2,((varw+S22)/B*(v-muv)-A/(B*B)*(v-muv)*F-(covvw+S12)/B*(w-muw)+C/(B*B)*(w-muw)*F-D/(B*B)*(x-mux)*F)*(v-muv)*A/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(varw+S22)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*A/(B*B)*F-((varw+S22)/B*(v-muv)-A/(B*B)*(v-muv)*F-(covvw+S12)/B*(w-muw)+C/(B*B)*(w-muw)*F-D/(B*B)*(x-mux)*F)*(w-muw)*C/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(covvw+S12)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*C/(B*B)*F+((varw+S22)/B*(v-muv)-A/(B*B)*(v-muv)*F-(covvw+S12)/B*(w-muw)+C/(B*B)*(w-muw)*F-D/(B*B)*(x-mux)*F)*(x-mux)*D/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*D/(B*B)*F-(varw+S22)/B+A/(B*B)*F);
    j.set(0,3,(-A/(B*B)*(v-muv)*G-(varx+S33)/B*(w-muw)+C/(B*B)*(w-muw)*G+(covwx+S23)/B*(x-mux)-D/(B*B)*(x-mux)*G)*(v-muv)*A/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*A/(B*B)*G-(-A/(B*B)*(v-muv)*G-(varx+S33)/B*(w-muw)+C/(B*B)*(w-muw)*G+(covwx+S23)/B*(x-mux)-D/(B*B)*(x-mux)*G)*(w-muw)*C/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(varx+S33)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*C/(B*B)*G+(-A/(B*B)*(v-muv)*G-(varx+S33)/B*(w-muw)+C/(B*B)*(w-muw)*G+(covwx+S23)/B*(x-mux)-D/(B*B)*(x-mux)*G)*(x-mux)*D/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(covwx+S23)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*D/(B*B)*G+A/(B*B)*G);
    j.set(0,4,(-A/(B*B)*(v-muv)*H-(-covwx-S23)/B*(w-muw)+C/(B*B)*(w-muw)*H+(-varw-S22)/B*(x-mux)-D/(B*B)*(x-mux)*H)*(v-muv)*A/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*A/(B*B)*H-(-A/(B*B)*(v-muv)*H-(-covwx-S23)/B*(w-muw)+C/(B*B)*(w-muw)*H+(-varw-S22)/B*(x-mux)-D/(B*B)*(x-mux)*H)*(w-muw)*C/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(-covwx-S23)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*C/(B*B)*H+(-A/(B*B)*(v-muv)*H-(-covwx-S23)/B*(w-muw)+C/(B*B)*(w-muw)*H+(-varw-S22)/B*(x-mux)-D/(B*B)*(x-mux)*H)*(x-mux)*D/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(-varw-S22)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*D/(B*B)*H+A/(B*B)*H);
    j.set(0,5,((-2*covwx-2*S23)/B*(v-muv)-A/(B*B)*(v-muv)*I-(-covvx-S13)/B*(w-muw)+C/(B*B)*(w-muw)*I+(covvw+S12)/B*(x-mux)-D/(B*B)*(x-mux)*I)*(v-muv)*A/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(-2*covwx-2*S23)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*A/(B*B)*I-((-2*covwx-2*S23)/B*(v-muv)-A/(B*B)*(v-muv)*I-(-covvx-S13)/B*(w-muw)+C/(B*B)*(w-muw)*I+(covvw+S12)/B*(x-mux)-D/(B*B)*(x-mux)*I)*(w-muw)*C/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(-covvx-S13)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*C/(B*B)*I+((-2*covwx-2*S23)/B*(v-muv)-A/(B*B)*(v-muv)*I-(-covvx-S13)/B*(w-muw)+C/(B*B)*(w-muw)*I+(covvw+S12)/B*(x-mux)-D/(B*B)*(x-mux)*I)*(x-mux)*D/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(covvw+S12)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*D/(B*B)*I-(-2*covwx-2*S23)/B+A/(B*B)*I);
    j.set(1,0,-(C/(B*B)*(v-muv)*A+(varx+S33)/B*(w-muw)-E/(B*B)*(w-muw)*A-(covwx+S23)/B*(x-mux)+J/(B*B)*(x-mux)*A)*(v-muv)*C/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*C/(B*B)*A+(C/(B*B)*(v-muv)*A+(varx+S33)/B*(w-muw)-E/(B*B)*(w-muw)*A-(covwx+S23)/B*(x-mux)+J/(B*B)*(x-mux)*A)*(w-muw)*E/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*(varx+S33)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*E/(B*B)*A-(C/(B*B)*(v-muv)*A+(varx+S33)/B*(w-muw)-E/(B*B)*(w-muw)*A-(covwx+S23)/B*(x-mux)+J/(B*B)*(x-mux)*A)*(x-mux)*J/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*(covwx+S23)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*J/(B*B)*A-(varx+S33)/B+A/(B*B)*E);
    j.set(1,1,-(C/(B*B)*(v-muv)*E-(E*E)/(B*B)*(w-muw)+J/(B*B)*(x-mux)*E)*(v-muv)*C/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*C/(B*B)*E+(C/(B*B)*(v-muv)*E-(E*E)/(B*B)*(w-muw)+J/(B*B)*(x-mux)*E)*(w-muw)*E/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*(E*E)/(B*B)-(C/(B*B)*(v-muv)*E-(E*E)/(B*B)*(w-muw)+J/(B*B)*(x-mux)*E)*(x-mux)*J/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*J/(B*B)*E+(E*E)/(B*B));
    j.set(1,2,-(-(covvw+S12)/B*(v-muv)+C/(B*B)*(v-muv)*F+(varv+S11)/B*(w-muw)-E/(B*B)*(w-muw)*F+J/(B*B)*(x-mux)*F)*(v-muv)*C/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*(covvw+S12)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*C/(B*B)*F+(-(covvw+S12)/B*(v-muv)+C/(B*B)*(v-muv)*F+(varv+S11)/B*(w-muw)-E/(B*B)*(w-muw)*F+J/(B*B)*(x-mux)*F)*(w-muw)*E/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*(varv+S11)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*E/(B*B)*F-(-(covvw+S12)/B*(v-muv)+C/(B*B)*(v-muv)*F+(varv+S11)/B*(w-muw)-E/(B*B)*(w-muw)*F+J/(B*B)*(x-mux)*F)*(x-mux)*J/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*J/(B*B)*F-(varv+S11)/B+E/(B*B)*F);
    j.set(1,3,-(-(varx+S33)/B*(v-muv)+C/(B*B)*(v-muv)*G-E/(B*B)*(w-muw)*G-(-covvx-S13)/B*(x-mux)+J/(B*B)*(x-mux)*G)*(v-muv)*C/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*(varx+S33)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*C/(B*B)*G+(-(varx+S33)/B*(v-muv)+C/(B*B)*(v-muv)*G-E/(B*B)*(w-muw)*G-(-covvx-S13)/B*(x-mux)+J/(B*B)*(x-mux)*G)*(w-muw)*E/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*E/(B*B)*G-(-(varx+S33)/B*(v-muv)+C/(B*B)*(v-muv)*G-E/(B*B)*(w-muw)*G-(-covvx-S13)/B*(x-mux)+J/(B*B)*(x-mux)*G)*(x-mux)*J/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*(-covvx-S13)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*J/(B*B)*G+E/(B*B)*G);
    j.set(1,4,-(-(-covwx-S23)/B*(v-muv)+C/(B*B)*(v-muv)*H+(-2*covvx-2*S13)/B*(w-muw)-E/(B*B)*(w-muw)*H-(-covvw-S12)/B*(x-mux)+J/(B*B)*(x-mux)*H)*(v-muv)*C/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*(-covwx-S23)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*C/(B*B)*H+(-(-covwx-S23)/B*(v-muv)+C/(B*B)*(v-muv)*H+(-2*covvx-2*S13)/B*(w-muw)-E/(B*B)*(w-muw)*H-(-covvw-S12)/B*(x-mux)+J/(B*B)*(x-mux)*H)*(w-muw)*E/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*(-2*covvx-2*S13)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*E/(B*B)*H-(-(-covwx-S23)/B*(v-muv)+C/(B*B)*(v-muv)*H+(-2*covvx-2*S13)/B*(w-muw)-E/(B*B)*(w-muw)*H-(-covvw-S12)/B*(x-mux)+J/(B*B)*(x-mux)*H)*(x-mux)*J/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*(-covvw-S12)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*J/(B*B)*H-(-2*covvx-2*S13)/B+E/(B*B)*H);
    j.set(1,5,-(-(-covvx-S13)/B*(v-muv)+C/(B*B)*(v-muv)*I-E/(B*B)*(w-muw)*I-(varv+S11)/B*(x-mux)+J/(B*B)*(x-mux)*I)*(v-muv)*C/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*(-covvx-S13)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*C/(B*B)*I+(-(-covvx-S13)/B*(v-muv)+C/(B*B)*(v-muv)*I-E/(B*B)*(w-muw)*I-(varv+S11)/B*(x-mux)+J/(B*B)*(x-mux)*I)*(w-muw)*E/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*E/(B*B)*I-(-(-covvx-S13)/B*(v-muv)+C/(B*B)*(v-muv)*I-E/(B*B)*(w-muw)*I-(varv+S11)/B*(x-mux)+J/(B*B)*(x-mux)*I)*(x-mux)*J/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*(varv+S11)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*J/(B*B)*I+E/(B*B)*I);
    j.set(2,0,(-D/(B*B)*(v-muv)*A-(covwx+S23)/B*(w-muw)+J/(B*B)*(w-muw)*A+(varw+S22)/B*(x-mux)-F/(B*B)*(x-mux)*A)*(v-muv)*D/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*D/(B*B)*A-(-D/(B*B)*(v-muv)*A-(covwx+S23)/B*(w-muw)+J/(B*B)*(w-muw)*A+(varw+S22)/B*(x-mux)-F/(B*B)*(x-mux)*A)*(w-muw)*J/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*(covwx+S23)/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*J/(B*B)*A+(-D/(B*B)*(v-muv)*A-(covwx+S23)/B*(w-muw)+J/(B*B)*(w-muw)*A+(varw+S22)/B*(x-mux)-F/(B*B)*(x-mux)*A)*(x-mux)*F/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(x-mux)*(varw+S22)/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(x-mux)*F/(B*B)*A-(varw+S22)/B+A/(B*B)*F);
    j.set(2,1,((-covvx-S13)/B*(v-muv)-D/(B*B)*(v-muv)*E+J/(B*B)*(w-muw)*E+(varv+S11)/B*(x-mux)-F/(B*B)*(x-mux)*E)*(v-muv)*D/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*(-covvx-S13)/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*D/(B*B)*E-((-covvx-S13)/B*(v-muv)-D/(B*B)*(v-muv)*E+J/(B*B)*(w-muw)*E+(varv+S11)/B*(x-mux)-F/(B*B)*(x-mux)*E)*(w-muw)*J/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*J/(B*B)*E+((-covvx-S13)/B*(v-muv)-D/(B*B)*(v-muv)*E+J/(B*B)*(w-muw)*E+(varv+S11)/B*(x-mux)-F/(B*B)*(x-mux)*E)*(x-mux)*F/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(x-mux)*(varv+S11)/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(x-mux)*F/(B*B)*E-(varv+S11)/B+E/(B*B)*F);
    j.set(2,2,(-D/(B*B)*(v-muv)*F+J/(B*B)*(w-muw)*F-(F*F)/(B*B)*(x-mux))*(v-muv)*D/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*D/(B*B)*F-(-D/(B*B)*(v-muv)*F+J/(B*B)*(w-muw)*F-(F*F)/(B*B)*(x-mux))*(w-muw)*J/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*J/(B*B)*F+(-D/(B*B)*(v-muv)*F+J/(B*B)*(w-muw)*F-(F*F)/(B*B)*(x-mux))*(x-mux)*F/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(x-mux)*(F*F)/(B*B)+(F*F)/(B*B));
    j.set(2,3,((covwx+S23)/B*(v-muv)-D/(B*B)*(v-muv)*G-(-covvx-S13)/B*(w-muw)+J/(B*B)*(w-muw)*G+(-2*covvw-2*S12)/B*(x-mux)-F/(B*B)*(x-mux)*G)*(v-muv)*D/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*(covwx+S23)/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*D/(B*B)*G-((covwx+S23)/B*(v-muv)-D/(B*B)*(v-muv)*G-(-covvx-S13)/B*(w-muw)+J/(B*B)*(w-muw)*G+(-2*covvw-2*S12)/B*(x-mux)-F/(B*B)*(x-mux)*G)*(w-muw)*J/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*(-covvx-S13)/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*J/(B*B)*G+((covwx+S23)/B*(v-muv)-D/(B*B)*(v-muv)*G-(-covvx-S13)/B*(w-muw)+J/(B*B)*(w-muw)*G+(-2*covvw-2*S12)/B*(x-mux)-F/(B*B)*(x-mux)*G)*(x-mux)*F/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(x-mux)*(-2*covvw-2*S12)/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(x-mux)*F/(B*B)*G-(-2*covvw-2*S12)/B+F/(B*B)*G);
    j.set(2,4,((-varw-S22)/B*(v-muv)-D/(B*B)*(v-muv)*H-(-covvw-S12)/B*(w-muw)+J/(B*B)*(w-muw)*H-F/(B*B)*(x-mux)*H)*(v-muv)*D/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*(-varw-S22)/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*D/(B*B)*H-((-varw-S22)/B*(v-muv)-D/(B*B)*(v-muv)*H-(-covvw-S12)/B*(w-muw)+J/(B*B)*(w-muw)*H-F/(B*B)*(x-mux)*H)*(w-muw)*J/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*(-covvw-S12)/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*J/(B*B)*H+((-varw-S22)/B*(v-muv)-D/(B*B)*(v-muv)*H-(-covvw-S12)/B*(w-muw)+J/(B*B)*(w-muw)*H-F/(B*B)*(x-mux)*H)*(x-mux)*F/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(x-mux)*F/(B*B)*H+F/(B*B)*H);
    j.set(2,5,((covvw+S12)/B*(v-muv)-D/(B*B)*(v-muv)*I-(varv+S11)/B*(w-muw)+J/(B*B)*(w-muw)*I-F/(B*B)*(x-mux)*I)*(v-muv)*D/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*(covvw+S12)/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(v-muv)*D/(B*B)*I-((covvw+S12)/B*(v-muv)-D/(B*B)*(v-muv)*I-(varv+S11)/B*(w-muw)+J/(B*B)*(w-muw)*I-F/(B*B)*(x-mux)*I)*(w-muw)*J/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*(varv+S11)/B+(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(w-muw)*J/(B*B)*I+((covvw+S12)/B*(v-muv)-D/(B*B)*(v-muv)*I-(varv+S11)/B*(w-muw)+J/(B*B)*(w-muw)*I-F/(B*B)*(x-mux)*I)*(x-mux)*F/B-(D/B*(v-muv)-J/B*(w-muw)+F/B*(x-mux))*(x-mux)*F/(B*B)*I+F/(B*B)*I);
    j.set(3,0,-(-(A*A)/(B*B)*(v-muv)+C/(B*B)*(w-muw)*A-D/(B*B)*(x-mux)*A)*(v-muv)*C/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*C/(B*B)*A+(-(A*A)/(B*B)*(v-muv)+C/(B*B)*(w-muw)*A-D/(B*B)*(x-mux)*A)*(w-muw)*E/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(varx+S33)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*E/(B*B)*A-(-(A*A)/(B*B)*(v-muv)+C/(B*B)*(w-muw)*A-D/(B*B)*(x-mux)*A)*(x-mux)*J/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(covwx+S23)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*J/(B*B)*A-C/(B*B)*A);
    j.set(3,1,-((varx+S33)/B*(v-muv)-A/(B*B)*(v-muv)*E+C/(B*B)*(w-muw)*E+(-covvx-S13)/B*(x-mux)-D/(B*B)*(x-mux)*E)*(v-muv)*C/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*C/(B*B)*E+((varx+S33)/B*(v-muv)-A/(B*B)*(v-muv)*E+C/(B*B)*(w-muw)*E+(-covvx-S13)/B*(x-mux)-D/(B*B)*(x-mux)*E)*(w-muw)*E/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(E*E)/(B*B)-((varx+S33)/B*(v-muv)-A/(B*B)*(v-muv)*E+C/(B*B)*(w-muw)*E+(-covvx-S13)/B*(x-mux)-D/(B*B)*(x-mux)*E)*(x-mux)*J/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*J/(B*B)*E-C/(B*B)*E);
    j.set(3,2,-((varw+S22)/B*(v-muv)-A/(B*B)*(v-muv)*F-(covvw+S12)/B*(w-muw)+C/(B*B)*(w-muw)*F-D/(B*B)*(x-mux)*F)*(v-muv)*C/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(covvw+S12)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*C/(B*B)*F+((varw+S22)/B*(v-muv)-A/(B*B)*(v-muv)*F-(covvw+S12)/B*(w-muw)+C/(B*B)*(w-muw)*F-D/(B*B)*(x-mux)*F)*(w-muw)*E/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(varv+S11)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*E/(B*B)*F-((varw+S22)/B*(v-muv)-A/(B*B)*(v-muv)*F-(covvw+S12)/B*(w-muw)+C/(B*B)*(w-muw)*F-D/(B*B)*(x-mux)*F)*(x-mux)*J/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*J/(B*B)*F+(covvw+S12)/B-C/(B*B)*F);
    j.set(3,3,-(-A/(B*B)*(v-muv)*G-(varx+S33)/B*(w-muw)+C/(B*B)*(w-muw)*G+(covwx+S23)/B*(x-mux)-D/(B*B)*(x-mux)*G)*(v-muv)*C/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(varx+S33)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*C/(B*B)*G+(-A/(B*B)*(v-muv)*G-(varx+S33)/B*(w-muw)+C/(B*B)*(w-muw)*G+(covwx+S23)/B*(x-mux)-D/(B*B)*(x-mux)*G)*(w-muw)*E/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*E/(B*B)*G-(-A/(B*B)*(v-muv)*G-(varx+S33)/B*(w-muw)+C/(B*B)*(w-muw)*G+(covwx+S23)/B*(x-mux)-D/(B*B)*(x-mux)*G)*(x-mux)*J/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(-covvx-S13)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*J/(B*B)*G+(varx+S33)/B-C/(B*B)*G);
    j.set(3,4,-(-A/(B*B)*(v-muv)*H-(-covwx-S23)/B*(w-muw)+C/(B*B)*(w-muw)*H+(-varw-S22)/B*(x-mux)-D/(B*B)*(x-mux)*H)*(v-muv)*C/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(-covwx-S23)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*C/(B*B)*H+(-A/(B*B)*(v-muv)*H-(-covwx-S23)/B*(w-muw)+C/(B*B)*(w-muw)*H+(-varw-S22)/B*(x-mux)-D/(B*B)*(x-mux)*H)*(w-muw)*E/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(-2*covvx-2*S13)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*E/(B*B)*H-(-A/(B*B)*(v-muv)*H-(-covwx-S23)/B*(w-muw)+C/(B*B)*(w-muw)*H+(-varw-S22)/B*(x-mux)-D/(B*B)*(x-mux)*H)*(x-mux)*J/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(-covvw-S12)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*J/(B*B)*H+(-covwx-S23)/B-C/(B*B)*H);
    j.set(3,5,-((-2*covwx-2*S23)/B*(v-muv)-A/(B*B)*(v-muv)*I-(-covvx-S13)/B*(w-muw)+C/(B*B)*(w-muw)*I+(covvw+S12)/B*(x-mux)-D/(B*B)*(x-mux)*I)*(v-muv)*C/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(-covvx-S13)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*C/(B*B)*I+((-2*covwx-2*S23)/B*(v-muv)-A/(B*B)*(v-muv)*I-(-covvx-S13)/B*(w-muw)+C/(B*B)*(w-muw)*I+(covvw+S12)/B*(x-mux)-D/(B*B)*(x-mux)*I)*(w-muw)*E/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*E/(B*B)*I-((-2*covwx-2*S23)/B*(v-muv)-A/(B*B)*(v-muv)*I-(-covvx-S13)/B*(w-muw)+C/(B*B)*(w-muw)*I+(covvw+S12)/B*(x-mux)-D/(B*B)*(x-mux)*I)*(x-mux)*J/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(varv+S11)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*J/(B*B)*I+(-covvx-S13)/B-C/(B*B)*I);
    j.set(4,0,(-(A*A)/(B*B)*(v-muv)+C/(B*B)*(w-muw)*A-D/(B*B)*(x-mux)*A)*(v-muv)*D/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*D/(B*B)*A-(-(A*A)/(B*B)*(v-muv)+C/(B*B)*(w-muw)*A-D/(B*B)*(x-mux)*A)*(w-muw)*J/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(covwx+S23)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*J/(B*B)*A+(-(A*A)/(B*B)*(v-muv)+C/(B*B)*(w-muw)*A-D/(B*B)*(x-mux)*A)*(x-mux)*F/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(varw+S22)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*F/(B*B)*A+D/(B*B)*A);
    j.set(4,1,((varx+S33)/B*(v-muv)-A/(B*B)*(v-muv)*E+C/(B*B)*(w-muw)*E+(-covvx-S13)/B*(x-mux)-D/(B*B)*(x-mux)*E)*(v-muv)*D/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(-covvx-S13)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*D/(B*B)*E-((varx+S33)/B*(v-muv)-A/(B*B)*(v-muv)*E+C/(B*B)*(w-muw)*E+(-covvx-S13)/B*(x-mux)-D/(B*B)*(x-mux)*E)*(w-muw)*J/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*J/(B*B)*E+((varx+S33)/B*(v-muv)-A/(B*B)*(v-muv)*E+C/(B*B)*(w-muw)*E+(-covvx-S13)/B*(x-mux)-D/(B*B)*(x-mux)*E)*(x-mux)*F/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(varv+S11)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*F/(B*B)*E-(-covvx-S13)/B+D/(B*B)*E);
    j.set(4,2,((varw+S22)/B*(v-muv)-A/(B*B)*(v-muv)*F-(covvw+S12)/B*(w-muw)+C/(B*B)*(w-muw)*F-D/(B*B)*(x-mux)*F)*(v-muv)*D/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*D/(B*B)*F-((varw+S22)/B*(v-muv)-A/(B*B)*(v-muv)*F-(covvw+S12)/B*(w-muw)+C/(B*B)*(w-muw)*F-D/(B*B)*(x-mux)*F)*(w-muw)*J/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*J/(B*B)*F+((varw+S22)/B*(v-muv)-A/(B*B)*(v-muv)*F-(covvw+S12)/B*(w-muw)+C/(B*B)*(w-muw)*F-D/(B*B)*(x-mux)*F)*(x-mux)*F/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(F*F)/(B*B)+D/(B*B)*F);
    j.set(4,3,(-A/(B*B)*(v-muv)*G-(varx+S33)/B*(w-muw)+C/(B*B)*(w-muw)*G+(covwx+S23)/B*(x-mux)-D/(B*B)*(x-mux)*G)*(v-muv)*D/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(covwx+S23)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*D/(B*B)*G-(-A/(B*B)*(v-muv)*G-(varx+S33)/B*(w-muw)+C/(B*B)*(w-muw)*G+(covwx+S23)/B*(x-mux)-D/(B*B)*(x-mux)*G)*(w-muw)*J/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(-covvx-S13)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*J/(B*B)*G+(-A/(B*B)*(v-muv)*G-(varx+S33)/B*(w-muw)+C/(B*B)*(w-muw)*G+(covwx+S23)/B*(x-mux)-D/(B*B)*(x-mux)*G)*(x-mux)*F/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*(-2*covvw-2*S12)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*F/(B*B)*G-(covwx+S23)/B+D/(B*B)*G);
    j.set(4,4,(-A/(B*B)*(v-muv)*H-(-covwx-S23)/B*(w-muw)+C/(B*B)*(w-muw)*H+(-varw-S22)/B*(x-mux)-D/(B*B)*(x-mux)*H)*(v-muv)*D/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(-varw-S22)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*D/(B*B)*H-(-A/(B*B)*(v-muv)*H-(-covwx-S23)/B*(w-muw)+C/(B*B)*(w-muw)*H+(-varw-S22)/B*(x-mux)-D/(B*B)*(x-mux)*H)*(w-muw)*J/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(-covvw-S12)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*J/(B*B)*H+(-A/(B*B)*(v-muv)*H-(-covwx-S23)/B*(w-muw)+C/(B*B)*(w-muw)*H+(-varw-S22)/B*(x-mux)-D/(B*B)*(x-mux)*H)*(x-mux)*F/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*F/(B*B)*H-(-varw-S22)/B+D/(B*B)*H);
    j.set(4,5,((-2*covwx-2*S23)/B*(v-muv)-A/(B*B)*(v-muv)*I-(-covvx-S13)/B*(w-muw)+C/(B*B)*(w-muw)*I+(covvw+S12)/B*(x-mux)-D/(B*B)*(x-mux)*I)*(v-muv)*D/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*(covvw+S12)/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(v-muv)*D/(B*B)*I-((-2*covwx-2*S23)/B*(v-muv)-A/(B*B)*(v-muv)*I-(-covvx-S13)/B*(w-muw)+C/(B*B)*(w-muw)*I+(covvw+S12)/B*(x-mux)-D/(B*B)*(x-mux)*I)*(w-muw)*J/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*(varv+S11)/B+(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(w-muw)*J/(B*B)*I+((-2*covwx-2*S23)/B*(v-muv)-A/(B*B)*(v-muv)*I-(-covvx-S13)/B*(w-muw)+C/(B*B)*(w-muw)*I+(covvw+S12)/B*(x-mux)-D/(B*B)*(x-mux)*I)*(x-mux)*F/B-(A/B*(v-muv)-C/B*(w-muw)+D/B*(x-mux))*(x-mux)*F/(B*B)*I-(covvw+S12)/B+D/(B*B)*I);
    j.set(5,0,(C/(B*B)*(v-muv)*A+(varx+S33)/B*(w-muw)-E/(B*B)*(w-muw)*A-(covwx+S23)/B*(x-mux)+J/(B*B)*(x-mux)*A)*(v-muv)*D/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*D/(B*B)*A-(C/(B*B)*(v-muv)*A+(varx+S33)/B*(w-muw)-E/(B*B)*(w-muw)*A-(covwx+S23)/B*(x-mux)+J/(B*B)*(x-mux)*A)*(w-muw)*J/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*(covwx+S23)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*J/(B*B)*A+(C/(B*B)*(v-muv)*A+(varx+S33)/B*(w-muw)-E/(B*B)*(w-muw)*A-(covwx+S23)/B*(x-mux)+J/(B*B)*(x-mux)*A)*(x-mux)*F/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*(varw+S22)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*F/(B*B)*A+(covwx+S23)/B-J/(B*B)*A);
    j.set(5,1,(C/(B*B)*(v-muv)*E-(E*E)/(B*B)*(w-muw)+J/(B*B)*(x-mux)*E)*(v-muv)*D/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*(-covvx-S13)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*D/(B*B)*E-(C/(B*B)*(v-muv)*E-(E*E)/(B*B)*(w-muw)+J/(B*B)*(x-mux)*E)*(w-muw)*J/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*J/(B*B)*E+(C/(B*B)*(v-muv)*E-(E*E)/(B*B)*(w-muw)+J/(B*B)*(x-mux)*E)*(x-mux)*F/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*(varv+S11)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*F/(B*B)*E-J/(B*B)*E);
    j.set(5,2,(-(covvw+S12)/B*(v-muv)+C/(B*B)*(v-muv)*F+(varv+S11)/B*(w-muw)-E/(B*B)*(w-muw)*F+J/(B*B)*(x-mux)*F)*(v-muv)*D/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*D/(B*B)*F-(-(covvw+S12)/B*(v-muv)+C/(B*B)*(v-muv)*F+(varv+S11)/B*(w-muw)-E/(B*B)*(w-muw)*F+J/(B*B)*(x-mux)*F)*(w-muw)*J/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*J/(B*B)*F+(-(covvw+S12)/B*(v-muv)+C/(B*B)*(v-muv)*F+(varv+S11)/B*(w-muw)-E/(B*B)*(w-muw)*F+J/(B*B)*(x-mux)*F)*(x-mux)*F/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*(F*F)/(B*B)-J/(B*B)*F);
    j.set(5,3,(-(varx+S33)/B*(v-muv)+C/(B*B)*(v-muv)*G-E/(B*B)*(w-muw)*G-(-covvx-S13)/B*(x-mux)+J/(B*B)*(x-mux)*G)*(v-muv)*D/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*(covwx+S23)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*D/(B*B)*G-(-(varx+S33)/B*(v-muv)+C/(B*B)*(v-muv)*G-E/(B*B)*(w-muw)*G-(-covvx-S13)/B*(x-mux)+J/(B*B)*(x-mux)*G)*(w-muw)*J/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*(-covvx-S13)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*J/(B*B)*G+(-(varx+S33)/B*(v-muv)+C/(B*B)*(v-muv)*G-E/(B*B)*(w-muw)*G-(-covvx-S13)/B*(x-mux)+J/(B*B)*(x-mux)*G)*(x-mux)*F/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*(-2*covvw-2*S12)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*F/(B*B)*G+(-covvx-S13)/B-J/(B*B)*G);
    j.set(5,4,(-(-covwx-S23)/B*(v-muv)+C/(B*B)*(v-muv)*H+(-2*covvx-2*S13)/B*(w-muw)-E/(B*B)*(w-muw)*H-(-covvw-S12)/B*(x-mux)+J/(B*B)*(x-mux)*H)*(v-muv)*D/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*(-varw-S22)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*D/(B*B)*H-(-(-covwx-S23)/B*(v-muv)+C/(B*B)*(v-muv)*H+(-2*covvx-2*S13)/B*(w-muw)-E/(B*B)*(w-muw)*H-(-covvw-S12)/B*(x-mux)+J/(B*B)*(x-mux)*H)*(w-muw)*J/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*(-covvw-S12)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*J/(B*B)*H+(-(-covwx-S23)/B*(v-muv)+C/(B*B)*(v-muv)*H+(-2*covvx-2*S13)/B*(w-muw)-E/(B*B)*(w-muw)*H-(-covvw-S12)/B*(x-mux)+J/(B*B)*(x-mux)*H)*(x-mux)*F/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*F/(B*B)*H+(-covvw-S12)/B-J/(B*B)*H);
    j.set(5,5,(-(-covvx-S13)/B*(v-muv)+C/(B*B)*(v-muv)*I-E/(B*B)*(w-muw)*I-(varv+S11)/B*(x-mux)+J/(B*B)*(x-mux)*I)*(v-muv)*D/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*(covvw+S12)/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(v-muv)*D/(B*B)*I-(-(-covvx-S13)/B*(v-muv)+C/(B*B)*(v-muv)*I-E/(B*B)*(w-muw)*I-(varv+S11)/B*(x-mux)+J/(B*B)*(x-mux)*I)*(w-muw)*J/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*(varv+S11)/B+(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(w-muw)*J/(B*B)*I+(-(-covvx-S13)/B*(v-muv)+C/(B*B)*(v-muv)*I-E/(B*B)*(w-muw)*I-(varv+S11)/B*(x-mux)+J/(B*B)*(x-mux)*I)*(x-mux)*F/B-(-C/B*(v-muv)+E/B*(w-muw)-J/B*(x-mux))*(x-mux)*F/(B*B)*I+(varv+S11)/B-J/(B*B)*I);
    return j;
}

protected static int n1 = 20, n2 = 20;

}