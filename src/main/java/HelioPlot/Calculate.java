package HelioPlot;

import Jama.Matrix;
import java.text.DecimalFormat;
import java.util.TreeMap;

public class Calculate {

    static double age(double U, double Th, double Sm, double He) throws Error {
        double[] ae = Calculate.ageErr(U, 0, Th, 0, Sm, 0, He, 0);
        return ae[0];
    }

    /*
     * input: composition = U, sU, Th, sTh, Sm, sSm, He, sHe
     * output: ageErr and 1-sigma uncertainty, both in Ma
     */
    static double[] ageErr(double U, double sU, double Th, double sTh,
                        double Sm, double sSm, double He, double sHe) throws Error {
        double a, b, c, d, P, Lwm, dtdHe, dtdU, dtdTh, dtdSm, 
               dfdt, dfdHe, dfdU, dfdTh, dfdSm;
        double[] ae = new double[2];        
        if (Sm==Data.NAN){Sm = 0; sSm = 0;}
        if (U >=0 && Th >=0 && He >= 0 && Sm >= 0){
            a = 8*(137.88/138.88)*L238;
            b = (7/138.88)*L235;
            c = 6*L232;
            d = 0.1499*L147;
            P = (a+b)*U + c*Th + d*Sm;
            Lwm = ((a*L238 + b*L235)*U + c*L232*Th + d*L147*Sm)/P;
            ae[0] = (1/Lwm)*Math.log(1+(Lwm/P)*He);
            dfdt = a*Math.exp(L238*ae[0])*U + b*Math.exp(L235*ae[0])*U + 
                   c*Math.exp(L232*ae[0])*Th + d*Math.exp(L147*ae[0])*Sm;
            dfdHe = -1;
            dfdU = 8*(137.88/138.88)*(Math.exp(L238*ae[0])-1) +
                   (7/138.88)*(Math.exp(L235*ae[0])-1);
            dfdTh = 6*(Math.exp(L232*ae[0])-1);
            dfdSm = 0.1499*(Math.exp(L147*ae[0])-1);
            dtdHe = -dfdHe/dfdt;
            dtdU = -dfdU/dfdt;
            dtdTh = -dfdTh/dfdt;
            dtdSm = -dfdSm/dfdt;
            ae[1] = Math.sqrt(dtdU*dtdU*sU*sU + dtdTh*dtdTh*sTh*sTh + dtdSm*dtdSm*sSm*sSm + dtdHe*dtdHe*sHe*sHe);
        }
        return ae;
    }

    /*
     * input: composition = U, sU, Th, sTh, Sm, sSm, He, sHe, C
     * output: ageErr and 1-sigma uncertainty, both in Ma
     */
    static double[] ageErr(double[] composition) throws Error {
        double U, sU, Th, sTh, Sm, sSm, He, sHe;
        U = composition[0]; sU = composition[1];
        Th = composition[2]; sTh = composition[3];
        Sm = composition[4]; sSm = composition[5];
        He = composition[6]; sHe = composition[7];
        return ageErr(U, sU, Th, sTh, Sm, sSm, He, sHe);
    }

    static double[] ageErr(Sample sample) throws Exception{
        return ageErr(sample.U(),sample.sU(),sample.Th(),sample.sTh(),
                      sample.Sm(),sample.sSm(),sample.He(),sample.sHe());
    }

    static double[] logAgeErr(Sample sample) throws Exception{
        double[] ae = ageErr(sample),
                 logae = {Math.log(ae[0]),ae[1]/ae[0]};
        return logae;
    }

    /* input: array of two matrixes containing the logratio means and (co)variances
     * output: array with the age and standard deviation
    */
    static double[] ageErr(Matrix[] muCov) {
        double[] ae = new double[2];
        boolean doSm = muCov[0].getRowDimension() == 3;
        double dtdV, dtdW, dtdX=Data.NAN, V, W, X=Data.NAN, UHe, ThHe, SmHe;
        Matrix J = muCov[0].copy();
        V = muCov[0].get(0,0);
        W = muCov[0].get(1,0);
        if (doSm) X = muCov[0].get(2,0);
        UHe = Math.pow(10,V);
        ThHe = Math.pow(10,W);
        SmHe = doSm ? Math.pow(10,X) : 0d;
        ae[0] = Calculate.age(UHe, ThHe, SmHe, 1d);
        dtdV = Math.log(10)*(8*(1-Math.exp(L238*ae[0]))+7*(1-Math.exp(L235*ae[0]))/137.88)/
               (8*L238*Math.exp(L238*ae[0])+7*L235*Math.exp(L235*ae[0])/137.88);        
        dtdW = Math.log(10)*(1-Math.exp(L232*ae[0]))/(L232*Math.exp(L232*ae[0]));
        if (doSm) dtdX = Math.log(10)*(1-Math.exp(L147*ae[0]))/(L147*Math.exp(L147*ae[0]));
        J.set(0,0,dtdV);
        J.set(1,0,dtdW);
        if (doSm) J.set(2,0,dtdX);
        ae[1] = Math.sqrt(J.transpose().times(muCov[1]).times(J).get(0,0));
        return ae;
    }

    static double USmHet2Th(double U, double Sm, double He, double t) throws Error {
        return (He - (8*(137.88/138.88)*(Math.exp(t*L238)-1) +
               (7/138.88)*(Math.exp(t*L235)-1))*U -
               0.1499*(Math.exp(t*L147)-1)*Sm)/
               (6*(Math.exp(t*L232)-1));
    }

    static double ThSmHet2U(double Th, double Sm, double He, double t) throws Error {
        return (He - 6*(Math.exp(t*L232)-1)*Th -
             0.1499*(Math.exp(t*L147)-1)*Sm)/
             (8*(137.88/138.88)*(Math.exp(t*L238)-1) +
            (7/138.88)*(Math.exp(t*L235)-1));
    }

    static double tThSm2U(double t, double Th, double Sm) throws Error {
        double a = (8*(137.88/138.88)*(Math.exp(t*L238)-1) +
            (7/138.88)*(Math.exp(t*L235)-1)),
               b = 6*(Math.exp(t*L232)-1),
               c = 0.1499*(Math.exp(t*L147)-1),
               U = (1 - (b+1)*Th - (c+1)*Sm)/(a+1);
        return U;
    }

    static double tUSm2Th(double t, double U, double Sm) throws Error {
        double a = (8*(137.88/138.88)*(Math.exp(t*L238)-1) +
            (7/138.88)*(Math.exp(t*L235)-1)),
               b = 6*(Math.exp(t*L232)-1),
               c = 0.1499*(Math.exp(t*L147)-1),
               Th = (1 - (a+1)*U - (c+1)*Sm)/(b+1);
        return Th;
    }

    static double tHeSm2U(double t, double He, double Sm) throws Error {
        double a = (8*(137.88/138.88)*(Math.exp(t*L238)-1) +
            (7/138.88)*(Math.exp(t*L235)-1)),
               b = 6*(Math.exp(t*L232)-1),
               c = 0.1499*(Math.exp(t*L147)-1),
               U = (He + (b-c)*Sm - b)/(a-b);
        return U;
    }

    /* input: U, Th, Sm, t
     * output: He
     */
    static double He(double U, double Th, double Sm, double t) throws Error{
        double He = (8*(137.88/138.88)*(Math.exp(t*L238)-1) +
            (7/138.88)*(Math.exp(t*L235)-1))*U +
             6*(Math.exp(t*L232)-1)*Th +
             0.1499*(Math.exp(t*L147)-1)*Sm;
        return He;
    }

    static double[] normalise(double U, double Th, double He) throws Error {
        double[] XYZ = {U/(U+Th+He), Th/(U+Th+He), He/(U+Th+He)};
        return XYZ;
    }

    static double[] vw2UThHe(double[] vw) throws Error {
        double[] UThHe = new double[3];
        UThHe[0] = Math.pow(10,vw[0])/(Math.pow(10,vw[0])+Math.pow(10,vw[1])+1);
        UThHe[1] = Math.pow(10,vw[1])/(Math.pow(10,vw[0])+Math.pow(10,vw[1])+1);
        UThHe[2] = 1/(Math.pow(10,vw[0])+Math.pow(10,vw[1])+1);
        return UThHe;
    }

    static double[] vwx2UThSmHe(double[] vwx) throws Error {
        double[] UThSmHe = new double[4];
        double denom = Math.pow(10,vwx[0]) + Math.pow(10,vwx[1]) + 1;
        if (vwx[2] != Data.NAN) {denom += Math.pow(10,vwx[2]);}
        UThSmHe[0] = Math.pow(10,vwx[0])/denom;
        UThSmHe[1] = Math.pow(10,vwx[1])/denom;
        UThSmHe[2] = (vwx[2] != Data.NAN) ? Math.pow(10,vwx[2])/denom : 0;
        UThSmHe[3] = 1/denom;
        return UThSmHe;
    }

    static double[] vwx2UThSmHe(Matrix vwx) throws Error {
        double[] UThSmHe = new double[4];
        if (vwx.getRowDimension()==2){
            double[] UThHe = vw2UThHe(vwx.getRowPackedCopy());
            UThSmHe[0] = UThHe[0];
            UThSmHe[1] = UThHe[1];
            UThSmHe[2] = 0;
            UThSmHe[3] = UThHe[2];
        } else {
            UThSmHe = vwx2UThSmHe(vwx.getRowPackedCopy());
        }
        return UThSmHe;
    }

    static double[] UThHe2vw(double[] xyz) throws Error {
        double vw[] = new double[2];
        vw[0] = Math.log10(xyz[0]/xyz[2]);
        vw[1] = Math.log10(xyz[1]/xyz[2]);
        return vw;
    }

    static double[] UThHe2vw(double x, double y, double z) throws Error {
        double vw[] = new double[2];
        vw[0] = Math.log10(x/z);
        vw[1] = Math.log10(y/z);
        return vw;
    }

    static void printTree(TreeMap<Double,Double> tree) throws Error {
        double t = tree.firstKey();
        while (t != tree.lastKey()) {
            System.out.println(t + " " + tree.get(t));
            t = tree.higherKey(t);
        }
    }

    /*
     * Rounds the ageErr t to the nearest ka, Ma or Ga
     */
    static public double roundAge(double t) throws Error {
        if (t<1){
            return Math.ceil(t*1000)/1000;
        } else if (t>1000) {
            return Math.ceil(t/1000)*1000;
        } else {
            return Math.ceil(t);
        }
    }

    static public String formatAgeLabel(double t) throws Error {
        return formatAgeLabel(t,"#");
    }

    static public String formatAgeLabel(double t, String format) throws Error {
        DecimalFormat fm = new DecimalFormat(format);
        if (t<1){
            return (fm.format(t*1000) + " ka");
        } else if (t>1000) {
            return (fm.format(t/1000) + " Ga");
        } else {
            return (fm.format(t) + " Ma");
        }
    }

    static double L238 = 1.55125e-4, L235 = 9.8485e-4, L232 = 0.49475e-4, L147=6.54e-6;

}