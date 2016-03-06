package HelioPlot;

public interface Coloured {

    public void setMinRGB(int r, int g, int b);
    
    public int[] getMinRGB();
    
    public void setMaxRGB(int r, int g, int b);

    public int[] getMaxRGB();
    
    public int getNumColours();
    
    public int[][] getColours();
    
}
