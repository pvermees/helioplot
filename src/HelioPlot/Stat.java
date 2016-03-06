package HelioPlot;

import Jama.Matrix;
import java.util.ArrayList;

/*
*   Class   Stat
*
*   USAGE:  Statistical functions
*
*   WRITTEN BY: Dr Michael Thomas Flanagan
*
*   DATE:    June 2002 as part of Fmath
*   AMENDED: 12 May 2003 Statistics separated out from Fmath as a new class
*   DATE:    18 June 2005, 5 January 2006, 25 April 2006, 12, 21 November 2006
*            4 December 2006 (renaming of cfd and pdf methods - older version also retained)
*            31 December 2006, March 2007, 14 April 2007, 19 October 2007, 27 February 2008
*            29 march 2008, 7 April 2008, 29 April 2008 - 13 May 2008, 22-31 May 2008,
*            4-10 June 2008, 27 June 2008, 2-5 July 2008, 23 July 2008, 31 July 2008,
*            2-4 August 2008,  20 August 2008, 5-10 September 2008, 19 September 2008, 28 September 2008
*
*   DOCUMENTATION:
*   See Michael Thomas Flanagan's Java library on-line web page:
*   http://www.ee.ucl.ac.uk/~mflanaga/java/Stat.html
*   http://www.ee.ucl.ac.uk/~mflanaga/java/
*
*   Copyright (c) 2002 - 2008 Michael Thomas Flanagan
*
*   PERMISSION TO COPY:
*
* Permission to use, copy and modify this software and its documentation for NON-COMMERCIAL purposes is granted, without fee,
* provided that an acknowledgement to the author, Dr Michael Thomas Flanagan at www.ee.ucl.ac.uk/~mflanaga, appears in all copies
* and associated documentation or publications.
*
* Redistributions of the source code of this source code, or parts of the source codes, must retain the above copyright notice, this list of conditions
* and the following disclaimer and requires written permission from the Michael Thomas Flanagan:
*
* Redistribution in binary form of all or parts of this class must reproduce the above copyright notice, this list of conditions and
* the following disclaimer in the documentation and/or other materials provided with the distribution and requires written permission from the Michael Thomas Flanagan:
*
* Dr Michael Thomas Flanagan makes no representations about the suitability or fitness of the software for any or for a particular purpose.
* Dr Michael Thomas Flanagan shall not be liable for any damages suffered as a result of using, modifying or distributing this software
* or its derivatives.
*
***************************************************************************************/
import java.util.TreeMap;

public class Stat {

    protected ArrayList<Object> array = null;    // internal array

    // STATIC VARIABLES

    // A small number close to the smallest representable floating point number
    public static final double FPMIN = 1e-300;

    // PRIVATE MEMBERS FOR USE IN GAMMA FUNCTION METHODS AND HISTOGRAM CONSTRUCTION METHODS

    // GAMMA FUNCTIONS
    //  Lanczos Gamma Function approximation - N (number of coefficients -1)
    protected static int lgfN = 6;
    //  Lanczos Gamma Function approximation - Coefficients
    protected static double[] lgfCoeff = {1.000000000190015, 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179E-2, -0.5395239384953E-5};
    //  Lanczos Gamma Function approximation - small gamma
    protected static double lgfGamma = 5.0;
    //  Maximum number of iterations allowed in Incomplete Gamma Function calculations
    protected static int igfiter = 1000;
    //  Tolerance used in terminating series in Incomplete Gamma Function calculations
    protected static double igfeps = 1e-8;

    // CONSTRUCTORS
    public Stat(){
         this.array = new ArrayList<Object>();
    }

    // Chi-Square Cumulative Distribution Function
    // probability that an observed chi-square value for a correct model should be less than chiSquare
    // nu  =  the degrees of freedom
    public static double chiSquareCDF(double chiSquare, int nu){
            if(nu<=0)throw new IllegalArgumentException("The degrees of freedom [nu], " + nu + ", must be greater than zero");
            return Stat.incompleteGamma((double)nu/2.0D, chiSquare/2.0D);
    }

    // Regularised Incomplete Gamma Function P(a,x) = integral from zero to x of (exp(-t)t^(a-1))dt
    // Retained for backward compatibility
    public static double incompleteGamma(double a, double x){
        return regularisedGammaFunction(a, x);
    }

    // Regularised Incomplete Gamma Function P(a,x) = integral from zero to x of (exp(-t)t^(a-1))dt
    public static double regularizedGammaFunction(double a, double x){
        return regularisedGammaFunction(a, x);
    }

    // Regularised Incomplete Gamma Function P(a,x) = integral from zero to x of (exp(-t)t^(a-1))dt
    public static double regularisedGammaFunction(double a, double x){
            if(a<0.0D  || x<0.0D)throw new IllegalArgumentException("\nFunction defined only for a >= 0 and x>=0");
            double igf;

            if(x < a+1.0D){
                    // Series representation
                    igf = incompleteGammaSer(a, x);
            }
            else{
                    // Continued fraction representation
                    igf = incompleteGammaFract(a, x);
            }
            return igf;
    }      

    // Regularised Incomplete Gamma Function P(a,x) = integral from zero to x of (exp(-t)t^(a-1))dt
    // Series representation of the function - valid for x < a + 1
    public static double incompleteGammaSer(double a, double x){
            if(a<0.0D  || x<0.0D)throw new IllegalArgumentException("\nFunction defined only for a >= 0 and x>=0");
            if(x>=a+1) throw new IllegalArgumentException("\nx >= a+1   use Continued Fraction Representation");

            int i = 0;
            double igf = 0.0D;
            boolean check = true;

            double acopy = a;
            double sum = 1.0/a;
            double incr = sum;
            double loggamma = Stat.logGamma(a);

            while(check){
                    ++i;
                    ++a;
                    incr *= x/a;
                    sum += incr;
                    if(Math.abs(incr) < Math.abs(sum)*Stat.igfeps){
                            igf = sum*Math.exp(-x+acopy*Math.log(x)- loggamma);
                            check = false;
                    }
                    if(i>=Stat.igfiter){
                            check=false;
                            igf = sum*Math.exp(-x+acopy*Math.log(x)- loggamma);
                            System.out.println("\nMaximum number of iterations were exceeded in Stat.incompleteGammaSer().\nCurrent value returned.\nIncrement = "+String.valueOf(incr)+".\nSum = "+String.valueOf(sum)+".\nTolerance =  "+String.valueOf(igfeps));
                    }
            }
            return igf;
    }

    // Regularised Incomplete Gamma Function P(a,x) = integral from zero to x of (exp(-t)t^(a-1))dt
    // Continued Fraction representation of the function - valid for x >= a + 1
    // This method follows the general procedure used in Numerical Recipes for C,
    // The Art of Scientific Computing
    // by W H Press, S A Teukolsky, W T Vetterling & B P Flannery
    // Cambridge University Press,   http://www.nr.com/
    public static double incompleteGammaFract(double a, double x){
            if(a<0.0D  || x<0.0D)throw new IllegalArgumentException("\nFunction defined only for a >= 0 and x>=0");
            if(x<a+1) throw new IllegalArgumentException("\nx < a+1   Use Series Representation");

            int i = 0;
            double ii;
            double igf;
            boolean check = true;

            double loggamma = Stat.logGamma(a);
            double numer;
            double incr;
            double denom = x - a + 1.0D;
            double first = 1.0D/denom;
            double term = 1.0D/FPMIN;
            double prod = first;

            while(check){
                    ++i;
                    ii = (double)i;
                    numer = -ii*(ii - a);
                    denom += 2.0D;
                    first = numer*first + denom;
                    if(Math.abs(first) < Stat.FPMIN){
                        first = Stat.FPMIN;
                    }
                    term = denom + numer/term;
                    if(Math.abs(term) < Stat.FPMIN){
                        term = Stat.FPMIN;
                     }
                    first = 1.0D/first;
                    incr = first*term;
                    prod *= incr;
                    if(Math.abs(incr - 1.0D) < igfeps)check = false;
                    if(i>=Stat.igfiter){
                            check=false;
                            System.out.println("\nMaximum number of iterations were exceeded in Stat.incompleteGammaFract().\nCurrent value returned.\nIncrement - 1 = "+String.valueOf(incr-1)+".\nTolerance =  "+String.valueOf(igfeps));
                    }
            }
            igf = 1.0D - Math.exp(-x+a*Math.log(x)-loggamma)*prod;
            return igf;
    }

    // log to base e of the Gamma function
    // Lanczos approximation (6 terms)
    // Retained for backward compatibility
    public static double logGamma(double x){
            double xcopy = x;
            double fg;
            double first = x + lgfGamma + 0.5;
            double second = lgfCoeff[0];

            if(x>=0.0){
                    if(x>=1.0 && x-(int)x==0.0){
                            fg = Stat.logFactorial(x)-Math.log(x);
                    }
                    else{
                            first -= (x + 0.5)*Math.log(first);
                            for(int i=1; i<=lgfN; i++)second += lgfCoeff[i]/++xcopy;
                            fg = Math.log(Math.sqrt(2.0*Math.PI)*second/x) - first;
                    }
            }
            else{
                    fg = Math.PI/(Stat.gamma(1.0D-x)*Math.sin(Math.PI*x));

                    if(fg!=1.0/0.0 && fg!=-1.0/0.0){
                            if(fg<0){
                                     throw new IllegalArgumentException("\nThe gamma function is negative");
                            }
                            else{
                                    fg = Math.log(fg);
                            }
                    }
            }
            return fg;
    }

    // Gamma function
    // Lanczos approximation (6 terms)
    // retained for backward compatibity
    public static double gamma(double x){

            double xcopy = x;
            double first = x + lgfGamma + 0.5;
            double second = lgfCoeff[0];
            double fg;

            if(x>=0.0){
                    if(x>=1.0D && x-(int)x==0.0D){
                            fg = Stat.factorial(x)/x;
                    }
                    else{
                            first = Math.pow(first, x + 0.5)*Math.exp(-first);
                            for(int i=1; i<=lgfN; i++)second += lgfCoeff[i]/++xcopy;
                            fg = first*Math.sqrt(2.0*Math.PI)*second/x;
                    }
            }
            else{
                     fg = -Math.PI/(x*Stat.gamma(-x)*Math.sin(Math.PI*x));
            }
            return fg;
    }        

    // log to base e of the factorial of n
    // log[e](factorial) returned as double
    // numerical rounding may makes this an approximation
    public static double logFactorial(int n){
        if(n<0)throw new IllegalArgumentException("\nn, " + n + ", must be a positive integer\nIs a Gamma funtion [Fmath.gamma(x)] more appropriate?");
        double f = 0.0D;
        for(int i=2; i<=n; i++)f+=Math.log(i);
        return f;
    }

    // log to base e of the factorial of n
    // Argument is of type double but must be, numerically, an integer
    // log[e](factorial) returned as double
    // numerical rounding may makes this an approximation
    public static double logFactorial(long n){
        if(n<0)throw new IllegalArgumentException("\nn, " + n + ", must be a positive integer\nIs a Gamma funtion [Fmath.gamma(x)] more appropriate?");
        double f = 0.0D;
        long iCount = 2L;
        while(iCount<=n){
            f+=Math.log(iCount);
            iCount += 1L;
        }
        return f;
    }

    // log to base e of the factorial of n
    // Argument is of type double but must be, numerically, an integer
    // log[e](factorial) returned as double
    // numerical rounding may makes this an approximation
    public static double logFactorial(double n){
        if(n<0 || (n-Math.floor(n))!=0)throw new IllegalArgumentException("\nn must be a positive integer\nIs a Gamma funtion [Fmath.gamma(x)] more appropriate?");
        double f = 0.0D;
        double iCount = 2.0D;
        while(iCount<=n){
            f+=Math.log(iCount);
            iCount += 1.0D;
        }
        return f;
    }

    // factorial of n
    // Argument is of type double but must be, numerically, an integer
    // factorial returned as double but is, numerically, should be an integer
    // numerical rounding may makes this an approximation after n = 21
    public static double factorial(double n){
        if(n<0 || (n-Math.floor(n))!=0)throw new IllegalArgumentException("\nn must be a positive integer\nIs a Gamma funtion [Fmath.gamma(x)] more appropriate?");
        double f = 1.0D;
        double iCount = 2.0D;
        while(iCount<=n){
            f*=iCount;
            iCount += 1.0D;
        }
        return f;
    }        

     // Returns a binomial mass probabilty function
    public static double binomialPDF(double p, int n, int k){
            if(k<0 || n<0)throw new IllegalArgumentException("\nn and k must be greater than or equal to zero");
            if(k>n)throw new IllegalArgumentException("\nk is greater than n");
            return Math.floor(0.5D + Math.exp(Stat.logFactorial(n) - Stat.logFactorial(k) - Stat.logFactorial(n-k)))*Math.pow(p, k)*Math.pow(1.0D - p, n - k);
    }

    /* added by Pieter Vermeesch on 9 March 2009
     * double normalPDF(Matrix x, Matrix mu, Matrix cov)
     * x, mu = [nx1] vector
     * cov = [nxn] covariance matrix
     */
    public static double normalPDF(Matrix x, Matrix mu, Matrix cov){
        double out;
        double n = cov.getRowDimension(), m = cov.getColumnDimension();
        if (x.getRowDimension() != n | mu.getRowDimension() != n | n != m){return 0;}
        out = Math.exp(x.minus(mu).transpose().times(cov.inverse()).times(x.minus(mu)).times(-0.5).get(0, 0))
              / Math.sqrt(Math.pow(2*Math.PI,n)*cov.det());
        return out;
    }
        
    static public TreeMap<Double,Double> pdf2cdf(TreeMap<Double,Double> tree){
        double t = tree.firstKey(), newf, oldf;
        // cumulative sum
        while (t != tree.lastKey()) {
            oldf = tree.get(t);
            t = tree.higherKey(t);
            newf = tree.get(t);
            tree.put(t, oldf + newf);
        }
        // normalise
        t = tree.firstKey();
        double sum = tree.lastEntry().getValue();
        oldf = tree.get(t);
        tree.put(t, oldf/sum);
        while (t != tree.lastKey()) {
            t = tree.higherKey(t);
            oldf = tree.get(t);
            tree.put(t, oldf/sum);
        }
        return tree;
    }
    
    static double percentile(TreeMap<Double, Double> cdf, double p) throws Exception {
        double t = cdf.firstKey();
        while (cdf.get(t)<p) {
            t = cdf.higherKey(t);
        }
        return t;
    }

    static double printTree(TreeMap<Double, Double> f) throws Exception {
        double t = f.firstKey();
        while (f.get(t)!=f.lastKey()) {
            t = f.higherKey(t);
            System.out.println(t + " " + f.get(t));
        }
        return t;
    }

    static double mean(int n, Object[] x) throws Exception {
        double s = 0;
        for (int i=0; i<n; i++){
            s += (Double)x[i];
        }
        return s/n;
    }

    static double var(int n, Object[] x) throws Exception {
        double s2 = 0, m = mean(n,x);
        for (int i=0; i<n; i++){
            s2 += ((Double)x[i] - m)*((Double)x[i] - m);
        }
        return s2/(n-1);
    }

    static double cov(int n, Object[] x, Object[] y) throws Exception {
        double s2 = 0, mx = mean(n,x), my = mean(n,y);
        for (int i=0; i<n; i++){
            s2 += ((Double)x[i] - mx)*((Double)y[i] - mx);
        }
        return s2/(n-1);
    }
}