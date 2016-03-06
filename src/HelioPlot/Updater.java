/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package HelioPlot;

import java.awt.Desktop;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import javax.swing.JEditorPane;
import javax.swing.JOptionPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

/**
 *
 * @author pvermees
 */
public class Updater {
    
    public static void run() {
        try {
            String newversion = getLatestVersion();
            if (!newversion.equals(Main.VERSION)) showUpdateMessage();
        } catch (Exception e) {
            if (Main.DEBUGGINGMODE) {e.printStackTrace(System.out);}
        }
    }
    
    static String getURL() {
        return "http://helioplot.london-geochron.com";
    }

    static String getProgramName() {
        return "HelioPlot " + Main.VERSION;
    }    
    
    private static void showUpdateMessage() throws Exception {
        final String url = getURL(),
                     message = "<br>A new version is available. Would you like<br><br>" + 
                               "to download " + getProgramName() + 
                               " from <br><br><a href=\"\">" + url + "</a>?<br>";
        JEditorPane ep = myJEditorPane(message,url,true);
        int reply = JOptionPane.showOptionDialog(null, ep, "update available",  JOptionPane.YES_NO_OPTION, JOptionPane.INFORMATION_MESSAGE,null,null,null);        
        if (reply==0) openURL(url,true);
    }
    
    static JEditorPane myJEditorPane(String message, final String url, final boolean exit) {
        JEditorPane ep = new JEditorPane("text/html", "<html><body><center>" +
                            message + "</center></body></html>");
        // handle link events
        ep.addHyperlinkListener(new HyperlinkListener(){
        @Override
        public void hyperlinkUpdate(HyperlinkEvent e){
            if (e.getEventType().equals(HyperlinkEvent.EventType.ACTIVATED))
                openURL(url, exit);
            }
        });
        ep.setEditable(false);   
        return ep;
    }
    
    static void openURL(String url, boolean exit) {
        try {
            if(Desktop.isDesktopSupported()){
                Desktop.getDesktop().browse(new URI(url));
            } else {
                Runtime.getRuntime().exec("xdg-open " + url);
            }
            if (exit) System.exit(0);
        } catch (URISyntaxException e) {
            if (Main.DEBUGGINGMODE){e.printStackTrace(System.out);}
        } catch (IOException e) {
            if (Main.DEBUGGINGMODE){e.printStackTrace(System.out);}
        }
    }
    
    public static String getLatestVersion() throws Exception {
        String data = getData(versionURL),
               open = "[helioplot]",
               close = "[/helioplot]";
        int start = open.length();
        return data.substring(data.indexOf(open)+start,data.indexOf(close));
    }

    private static String getData(String address)throws Exception {
        URL url = new URL(address);        
        InputStream html;
        html = url.openStream();
        int c = 0;
        StringBuilder buffer = new StringBuilder("");
        while(c != -1) {
            c = html.read();            
            buffer.append((char)c);
        }
        return buffer.toString();
    }
    
    static final String versionURL = "http://ucl.ac.uk/~ucfbpve/software/version.html";
    
}

