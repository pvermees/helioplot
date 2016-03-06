package HelioPlot;

import java.awt.Color;
import java.util.ArrayList;

abstract public class ColouredObject implements Coloured {
    
    public ColouredObject(){
        try {
            this.colours = new int[3][this.numcolours];
            this.initColours();
        } catch (Exception ex) {
            if (Main.DEBUGGINGMODE){ex.printStackTrace(System.out);}
        }
    }
    
    protected ArrayList<Double> getAgeContours(double tmin, double tmax, int numcontours) throws Exception {
        ArrayList<Double> t = new ArrayList<Double>();
        t.add(Calculate.roundAge(tmin));
        for (int i=0; i<numcontours; i++){
            t.add(Calculate.roundAge(Math.exp(Math.log(tmin)+
                  i*(Math.log(tmax)-Math.log(tmin))/numcontours)));
        }
        t.add(Calculate.roundAge(tmax));
        return t;
    }
    
    // set the minimum end of the colour scale
    @Override
    public void setMinRGB(int r, int g, int b){
        rgb1[0] = r;
        rgb1[1] = g;
        rgb1[2] = b;
        this.initColours();
    }
    
    // get the minimum end of the colour scale (returns int[3])
    @Override
    public int[] getMinRGB(){
        return rgb1;
    }
    
    // set the maximum end of the colour scale
    @Override
    public void setMaxRGB(int r, int g, int b){
        rgb2[0] = r;
        rgb2[1] = g;
        rgb2[2] = b;
        this.initColours();
    }

    /*
     * get the minimum end of the colour scale (returns int[3])
     */
    @Override
    public int[] getMaxRGB(){
        return rgb2;
    }
    
    @Override
    public int getNumColours(){
        return numcolours;
    }
    
    @Override
    public int[][] getColours(){
        return colours;
    }
    
    /* 
     * returns an array of 3 integers between 0 and 255
     */
    public Color getColour(double val) throws Exception {
        int colournum;
        Color c;
        if (minC > 0 && maxC > 0){
            colournum = (int) Math.floor((numcolours-1)*
                          Math.abs((Math.log(val)-Math.log(minC))/
                                    (Math.log(maxC)-Math.log(minC))));
        } else {
            colournum = (int) Math.floor((numcolours-1)*
                          Math.abs((val-minC)/(maxC-minC)));
        }
        c = new Color(colours[0][colournum],colours[1][colournum],colours[2][colournum]);
        return c;
    }
    
    public final void initColours() {
        for (int i = 0; i<numcolours; i++){
            colours[0][i] = (int)(rgb1[0] + (rgb2[0] - rgb1[0])*i/numcolours);
            colours[1][i] = (int)(rgb1[1] + (rgb2[1] - rgb1[1])*i/numcolours);
            colours[2][i] = (int)(rgb1[2] + (rgb2[2] - rgb1[2])*i/numcolours);
        }
    }
    
    abstract void plot() throws Exception ;
    
    protected double minC, maxC;
    protected int numcolours = 100;
    protected int[][] colours;
    protected int[] rgb1 = {255,255,255}, rgb2 = {0,0,0};
    
}