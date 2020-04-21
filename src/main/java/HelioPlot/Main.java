package HelioPlot;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.HeadlessException;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import javax.swing.filechooser.FileFilter;
import java.text.DecimalFormat;
import javax.imageio.ImageIO;
import javax.swing.*;
import java.io.*;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumnModel;
import java.net.URL;
import java.util.ArrayList;
import javax.help.HelpBroker;
import javax.help.HelpSet;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.TableCellRenderer;
import org.freehep.graphicsbase.util.UserProperties;
import org.freehep.graphicsio.pdf.PDFGraphics2D;

public class Main extends JFrame implements TableModelListener {

    public Main() {
        super("Input");
        //this.setIconImage(new ImageIcon(getClass().getResource("icon.gif")).getImage());
        initComponents();
        data = new Data();
        // the following two lines are for debugging purposes only
        data.load("");
        this.loadData();
        helioplot = new HelioPlot(this);
        helioplot.setBounds(this.getBounds().x + this.getWidth(),
            this.getBounds().y, HelioPlot.SIZE, HelioPlot.SIZE);
        Updater.run();
    }

    private void initComponents(){
        setLookAndFeel();
        JButton button = new JButton("Plot!");
        button.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                plotbuttonEventHandler(e);
            }
        });
        this.add(button, BorderLayout.PAGE_END);
        this.pack();
        JMenuBar menubar = new JMenuBar();
        JMenu filemenu = new JMenu("File");
        JMenuItem newitem = new JMenuItem("New",KeyEvent.VK_N);
        newitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, ActionEvent.CTRL_MASK));
        newitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                newitemEventHandler(e);
            }
        });
        filemenu.add(newitem);
        JMenuItem openitem = new JMenuItem("Open",KeyEvent.VK_O);
        openitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, ActionEvent.CTRL_MASK));
        openitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                openitemEventHandler(e);
            }
        });
        filemenu.add(openitem);
        JMenuItem saveplotitem = new JMenuItem("Save Plot",KeyEvent.VK_S);
        saveplotitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, ActionEvent.CTRL_MASK));
        saveplotitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveplotEventHandler(e);
            }
        });
        filemenu.add(saveplotitem);
        JMenuItem savedataitem = new JMenuItem("Save Data",KeyEvent.VK_A);
        savedataitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A, ActionEvent.CTRL_MASK));
        savedataitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                savedataEventHandler(e);
            }
        });
        filemenu.add(savedataitem);
        JMenuItem exititem = new JMenuItem("Exit",KeyEvent.VK_W);
        exititem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_W, ActionEvent.CTRL_MASK));
        exititem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                exititemEventHandler(e);
            }
        });
        filemenu.add(exititem);
        menubar.add(filemenu);

        JMenu editmenu = new JMenu("Edit");
        JMenuItem copyitem = new JMenuItem("Copy",KeyEvent.VK_C);
        copyitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C, ActionEvent.CTRL_MASK));
        copyitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                copyitemEventHandler(e);
            }
        });
        editmenu.add(copyitem);
        JMenuItem pasteitem = new JMenuItem("Paste",KeyEvent.VK_V);
        pasteitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V, ActionEvent.CTRL_MASK));
        pasteitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                pasteitemEventHandler(e);
            }
        });
        editmenu.add(pasteitem);
        JMenuItem insertitem = new JMenuItem("Insert",KeyEvent.VK_I);
        insertitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_I, ActionEvent.CTRL_MASK));
        insertitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                insertitemEventHandler(e);
            }
        });
        editmenu.add(insertitem);
        JMenuItem deleteitem = new JMenuItem("Delete",KeyEvent.VK_D);
        deleteitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, ActionEvent.CTRL_MASK));
        deleteitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                deleteitemEventHandler(e);
            }
        });
        editmenu.add(deleteitem);
        menubar.add(editmenu);

        JMenu optionsmenu = new JMenu("Options");
        JMenu inputoptions = new JMenu("Input");
        ButtonGroup input = new ButtonGroup();
        uthhebutton = new JRadioButtonMenuItem("(U-Th)/He");
        input.add(uthhebutton);
        uthhebutton.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                uthhebuttonEventHandler(e);
            }
        });
        uthhebutton.setSelected(true);
        uthsmhebutton = new JRadioButtonMenuItem("(U-Th-Sm)/He");
        input.add(uthsmhebutton);
        uthsmhebutton.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                uthsmhebuttonEventHandler(e);
            }
        });
        inputoptions.add(uthhebutton);
        inputoptions.add(uthsmhebutton);
        optionsmenu.add(inputoptions);
        JMenu outputoptions = new JMenu("Output");
        ButtonGroup output = new ButtonGroup();
        ternarybutton = new JRadioButtonMenuItem("Ternary diagram");
        output.add(ternarybutton);
        ternarybutton.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // nothing to implement
            }
        });
        logratiobutton = new JRadioButtonMenuItem("Logratio plot");
        output.add(logratiobutton);
        logratiobutton.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // noting to implement
            }
        });
        outputoptions.add(ternarybutton);
        outputoptions.add(logratiobutton);
        optionsmenu.add(outputoptions);
        logratiobutton.setSelected(true);
        JMenu colouroptions = new JMenu("Colours");
        JMenuItem contourcolouritem = new JMenuItem("Contour Colour");
        contourcolouritem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                contourcolouritemEventHandler(e);
            }
        });
        colouroptions.add(contourcolouritem);
        JMenuItem ellipsecolouritem = new JMenuItem("Ellipse Colour");
        ellipsecolouritem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ellipsecolouritemEventHandler(e);
            }
        });
        colouroptions.add(ellipsecolouritem);
        optionsmenu.add(colouroptions);
        JMenuItem settingsitem = new JMenuItem("Settings");
        settingsitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                settingsitemEventHandler(e);
            }
        });
        optionsmenu.add(settingsitem);
        menubar.add(optionsmenu);

        JMenu helpmenu = new JMenu("Help");
        JMenuItem contentsitem = new JMenuItem("Contents",KeyEvent.VK_F1);
        contentsitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F1, ActionEvent.CTRL_MASK));
        contentsitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                contentsitemEventHandler(e);
            }
        });
        helpmenu.add(contentsitem);
        JMenuItem aboutitem = new JMenuItem("About");
        aboutitem.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                aboutitemEventHandler(e);
            }
        });
        helpmenu.add(aboutitem);
        menubar.add(helpmenu);

        this.setJMenuBar(menubar);

        samplename = new JTextField("sample name");
        samplename.setHorizontalAlignment(JTextField.CENTER);
        this.add(samplename,BorderLayout.PAGE_START);

        samplename.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                samplenameEnterHandler(evt);
            }
        });
        samplename.addFocusListener(new java.awt.event.FocusAdapter() {
            @Override
            public void focusLost(java.awt.event.FocusEvent evt) {
                samplenameFocusHandler(evt);
            }
        });

        model = new DefaultTableModel(100, numcols);
        model.addTableModelListener(this);
        /* this custom constructor changes the color of the last two columns
         */
        table = new JTable(model){
            @Override
            public Component prepareRenderer(TableCellRenderer renderer,
                                             int rowIndex, int vColIndex) {
                Component c = super.prepareRenderer(renderer, rowIndex, vColIndex);
                if (vColIndex >= numcols-2 && !isCellSelected(rowIndex, vColIndex)) {
                    c.setForeground(Color.blue);
                } else {
                    c.setForeground(Color.black);
                }
                return c;
            }
        };
        table.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
        table.setCellSelectionEnabled(true);
        scrollPane = new JScrollPane( table );
        lineTable = new LineNumberTable( table );
        clipboard = new ExcelAdapter(table, this);
        this.setTableModel(numcols);

        scrollPane.setRowHeaderView(lineTable);

        this.setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );
        this.getContentPane().add(scrollPane, BorderLayout.CENTER);

        table.putClientProperty( "JTable.autoStartsEdit", false);

        this.setSize(700, 500);
        this.setVisible(true);
    }

    public Data getData(){
        return data;
    }

    private void setTableModel(int n){
        numcols = n;
        int numrows = model.getRowCount();
        model = new DefaultTableModel(numrows, numcols){
            @Override public boolean isCellEditable(int row, int col) {
                return col < numcols-2;
            }
        };
        table.setModel(model);
        this.setColumnHeaders();
        model.addTableModelListener(this);
    }

    public boolean noColour(){
        try {
            ArrayList<Double> S = data.getSm(), 
                              C = data.getC();
            for (int i=0; i<C.size(); i++){
                if (S.get(i) != Data.NAN | C.get(i) != Data.NAN){return false;}
            }
            return true;
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
            return true;
        }
    }

    private void setColumnHeaders(){
        TableColumnModel colmod = lineTable.getMainTable().getColumnModel();
        colmod.getColumn(0).setHeaderValue("U");
        colmod.getColumn(1).setHeaderValue("σ(U)");
        colmod.getColumn(2).setHeaderValue("Th");
        colmod.getColumn(3).setHeaderValue("σ(Th)");
        if (numcols == N2){
            colmod.getColumn(4).setHeaderValue("Sm");
            colmod.getColumn(5).setHeaderValue("σ(Sm)");
            colmod.getColumn(6).setHeaderValue("He");
            colmod.getColumn(7).setHeaderValue("σ(He)");
            colmod.getColumn(8).setHeaderValue("t");
            colmod.getColumn(9).setHeaderValue("σ(t)");
        } else {
            colmod.getColumn(4).setHeaderValue("He");
            colmod.getColumn(5).setHeaderValue("σ(He)");
            colmod.getColumn(6).setHeaderValue("[C]");
            colmod.getColumn(7).setHeaderValue("t");
            colmod.getColumn(8).setHeaderValue("σ(t)");
        }
    }

    private void clearDisplay() {
        this.samplename.setText("");
        for (int i=0; i<numcols; i++){
            for (int j=0; j<table.getRowCount(); j++){
                table.setValueAt("", j, i);
            }
        }
    }

    /* sets the look to match the platform's OS (PC, Mac, Linux, ...)*/
    private static void setLookAndFeel(){
        try {
          UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch(ClassNotFoundException e) {
          throw new RuntimeException(e);
        } catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        } catch (InstantiationException e) {
            throw new RuntimeException(e);
        } catch (UnsupportedLookAndFeelException e) {
            throw new RuntimeException(e);
        }
    }

    private void uthhebuttonEventHandler(ActionEvent evt){
        try {
            setTableModel(N1);
            loadData();
            helioplot.getColourScale().setMinMaxC();
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    private void uthsmhebuttonEventHandler(ActionEvent evt){
        try {
            setTableModel(N2);
            loadData();
            helioplot.getColourScale().setMinMaxC();
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    private void plotbuttonEventHandler(ActionEvent evt) {
        try {
            if (pasted){
                helioplot.autoScale();
                pasted = false;
            }
            helioplot.setVisible(true);          
            helioplot.repaint();
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    private void newitemEventHandler(ActionEvent evt) {
        samplename.setText("");
        for (int i=0; i<table.getRowCount();i++){
            for (int j=0; j<table.getColumnCount();j++){
                table.setValueAt("", i, j);
            }
        }
    }

    private void openitemEventHandler(ActionEvent evt) {
        try {
            String filename;
            JFileChooser chooser;
            int returnVal;
            chooser = new JFileChooser(idir);
            returnVal = chooser.showOpenDialog((Component) evt.getSource());
            if ( returnVal == JFileChooser.APPROVE_OPTION ) {
                filename = chooser.getSelectedFile().getAbsolutePath();
                idir = chooser.getCurrentDirectory().getAbsolutePath();
                this.clearDisplay();
                data.load(filename);
                if (data.hasSm() && !data.hasColour()){
                    this.uthsmhebutton.doClick();
                }
                helioplot.autoScale();
                this.loadData();
            }
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    /* displays the composition in table and updates the plot */
    private void loadData() {
        try {
            this.samplename.setText(data.samplename);
            if (!data.samplename.equals("")){
                samplename.setText(data.samplename);
                for (int i = 0; i < data.size(); i++){
                    if (i>=model.getRowCount()){
                        model.addRow(new Object[]{""});
                    }
                    int[] cols = {0,1,2,3,4,5,6,7};
                    if (!this.UThSmHe()){
                        cols[4] = 6;
                        cols[5] = 7;
                        cols[6] = 8;
                    }
                    for (int j=0; j<numcols-2; j++){
                        if (data.get(i, cols[j]) == Data.NAN) {
                            model.setValueAt("", i, j);
                        } else {
                            model.setValueAt(data.get(i, cols[j]), i, j);
                        }
                    }
                }
            }
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    private void samplenameFocusHandler(java.awt.event.FocusEvent evt) {
        try {
            data.setSampleName(((JTextField)evt.getSource()).getText());
        } catch (Exception e){
            if (DEBUGGINGMODE){e.printStackTrace(System.out);}
        }
    }

    private void samplenameEnterHandler(java.awt.event.ActionEvent evt) {
        try {
            data.setSampleName(((JTextField)evt.getSource()).getText());
        } catch (Exception e){
            if (DEBUGGINGMODE){e.printStackTrace(System.out);}
        }
    }

    private void saveplotEventHandler(ActionEvent evt) {
        try {
            if (odir.equals("")) {
                odir = idir.equals("") ? System.getProperty("user.dir") : idir;
            }
            final JFileChooser fc = new JFileChooser();
            fc.setAcceptAllFileFilterUsed(false);
            fc.setCurrentDirectory(new File(odir));
            String png = "*.png, *.PNG", pdf = "*.pdf, *.PDF";
            FileNameExtensionFilter f1 = new FileNameExtensionFilter(pdf, "pdf", "PDF"),
                                    f2 = new FileNameExtensionFilter(png, "png", "PNG");           
            fc.addChoosableFileFilter(f1);
            fc.addChoosableFileFilter(f2);
            int returnVal = fc.showSaveDialog(this);
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                FileFilter fF = fc.getFileFilter();
                odir = fc.getSelectedFile().getParent() + File.separator;
                String fname = removeExtension(fc.getSelectedFile().getPath());
                if (fF.getDescription().equals(pdf)){
                    this.savePlot(odir+fname,".pdf");
                } else if (fF.getDescription().equals(png)) {
                    this.savePlot(odir+fname,".png");
                }
            } 
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    public void savePlot(String out, String extension) throws Exception {
        Dimension olddim = helioplot.plotpanel.getSize();
        helioplot.plotpanel.setSize(2*HelioPlot.SIZE,2*HelioPlot.SIZE);
        helioplot.plotpanel.revalidate();
        if (extension.equals(".pdf")) {
            PDFGraphics2D g = new PDFGraphics2D(new File(out+extension), helioplot.plotpanel.getSize());
            UserProperties p = new UserProperties();
            p.setProperty(PDFGraphics2D.TEXT_AS_SHAPES, false);
            g.setProperties(p);
            g.startExport();
            helioplot.plotpanel.print(g);
            g.endExport();
        } else if (extension.equals(".png")) {
            BufferedImage bi = new BufferedImage(helioplot.getPanelWidth(), 
                                                 helioplot.getPanelHeight(), 
                                                 BufferedImage.TYPE_INT_ARGB);
            Graphics2D g = bi.createGraphics();
            helioplot.plotpanel.print(g);
            ImageIO.write(bi, "PNG", new File(out+extension));            
        }
        helioplot.plotpanel.setSize(olddim);
        helioplot.plotpanel.revalidate();        
    }
    
    public static String removeExtension(String s) {

        String separator = System.getProperty("file.separator");
        String filename;

        // Remove the path upto the filename.
        int lastSeparatorIndex = s.lastIndexOf(separator);
        if (lastSeparatorIndex == -1) {
            filename = s;
        } else {
            filename = s.substring(lastSeparatorIndex + 1);
        }

        // Remove the extension.
        int extensionIndex = filename.lastIndexOf(".");
        if (extensionIndex == -1)
            return filename;

        return filename.substring(0, extensionIndex);
    }  

    private void savedataEventHandler(ActionEvent evt) {
        try {
            String initialDirectory = "";
            JFileChooser chooser;
            int returnVal;
            chooser = new JFileChooser(initialDirectory);
            FileFilter filter = new FileFilter() {
                @Override
                public boolean accept(File f) {
                    return f.getName().toLowerCase().endsWith(".csv") || f.isDirectory();
                }

                @Override
                public String getDescription() {
                    return "*.csv";
                }
            };
            chooser.setFileFilter(filter);
            chooser.setSelectedFile(new File("*.csv"));
            returnVal = chooser.showSaveDialog((Component) evt.getSource());
            if ( returnVal == JFileChooser.APPROVE_OPTION ) {
                String filepath = chooser.getSelectedFile().getAbsolutePath();
                data.write(filepath);
            }
        } catch (HeadlessException ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    private void copyitemEventHandler(ActionEvent evt) {
        clipboard.actionPerformed(evt);
    }

    private void pasteitemEventHandler(ActionEvent evt) {
        clipboard.actionPerformed(evt);
    }

    private void insertitemEventHandler(ActionEvent evt) {
        try {
            int r = table.getSelectedRow();
            if (r>-1){
                model.insertRow(r, new Object[]{""});
                lineTable.setModel(model);
                data.insertRow(r);
            }
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    private void deleteitemEventHandler(ActionEvent evt) {
        try {
            int r = table.getSelectedRow();
            if (r>-1){
                model.removeRow(r);
                data.removeRow(r);
            }
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    private void contourcolouritemEventHandler(ActionEvent evt) {
        try {
            if (ternarybutton.isSelected()){
                ColourChooser.createAndShowGUI(helioplot.getTernaryPlot());
            } else if (logratiobutton.isSelected()){
                ColourChooser.createAndShowGUI(helioplot.getLogRatioPlot());
            }
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    private void ellipsecolouritemEventHandler(ActionEvent evt) {
        ColourChooser.createAndShowGUI(helioplot.getColourScale());
    }

    private void settingsitemEventHandler(ActionEvent e) {
        Settings.createAndShowGUI(this);
    }

    @SuppressWarnings({"UseSpecificCatch", "BroadCatchBlock", "TooBroadCatch"})
    private void contentsitemEventHandler(ActionEvent evt) {
            // Identify the location of the .hs file
            String pathToHS = "/HelioHelp/docs/heliohelp-hs.xml";
        try {
            //Create a URL for the location of the help set
            hsURL = this.getClass().getResource(pathToHS);
            hs = new HelpSet(null, hsURL);
            // Create a HelpBroker object for manipulating the help set
            hb = hs.createHelpBroker();
            //Display help set
            hb.setDisplayed(true);
        } catch (Exception ee) {
            // Print info to the console if there is an exception
            System.out.println( "HelpSet " + ee.getMessage());
            System.out.println("Help Set "+ pathToHS +" not found");
        }
    }

    private void aboutitemEventHandler(ActionEvent evt) {
        String url = Updater.getURL();
        String message = Updater.getProgramName() + 
                    "<br><br>Pieter Vermeesch<br>" +
                    "<br>London Geochronology Centre<br>" + 
                    "<br><a href=\"\">" + url + "</a><br>";
        JEditorPane ep = Updater.myJEditorPane(message,url,false);
        JOptionPane.showMessageDialog(this, ep); 
    }

    private void exititemEventHandler(ActionEvent evt) {
        System.exit(0);
    }

    @Override
    public void tableChanged(TableModelEvent e) {
        DecimalFormat f = new DecimalFormat("#0.00");
        int c = e.getColumn();
        int r = e.getFirstRow();
        double val = 0.0;
        try {
            if (c<0) {
                // do nothing
            } else if (table.getValueAt(r, c).getClass().toString().equals("class java.lang.String")){
                val = Double.parseDouble((String) table.getValueAt(r, c));
            } else if ((table.getValueAt(r, c).getClass().toString().equals("class java.lang.Integer"))){
                val = (double) ((int)((Integer) table.getValueAt(r, c)));
            } else if ((table.getValueAt(r, c).getClass().toString().equals("class java.lang.Double"))){
                val = (double) ((Double) table.getValueAt(r, c));
            }
        } catch (NumberFormatException e1){
            val = Data.NAN;
        }
        try {
            // only column 8 (colour) is allowed to have zero or negative values
            if (this.UThSmHe() || c<6) {
                val = val > 0.0 ? val : Data.NAN;
            }
            // don't update the data when the last two columns are changed
            if (c>-1 & c<numcols-2){
                if (this.UThSmHe()) {
                    data.set(val,c,r);
                } else {
                    switch (c){
                        case 4: data.set(val,6,r); break;
                        case 5: data.set(val,7,r); break;
                        case 6: data.set(val,8,r); break;
                        default: data.set(val,c,r); break;
                    }
                }
                double[] AE = Calculate.ageErr(data.getRow(r));
                String age = "", err = "";
                if (AE[0] != Data.NAN){age = f.format(AE[0]);}
                if (AE[1] != Data.NAN){err = f.format(AE[1]);}
                model.setValueAt(age, r, numcols-2);
                model.setValueAt(err, r, numcols-1);
            }
        } catch (Exception ex){
            if (DEBUGGINGMODE) {ex.printStackTrace(System.out);}
        }
    }

    public boolean ternary(){
        return this.ternarybutton.isSelected();
    }

    public boolean UThSmHe(){
        return this.uthsmhebutton.isSelected();
    }

    public HelioPlot getHelioPlot(){
        return helioplot;
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                Main main = new Main();
            }
        });
    }

    private DefaultTableModel model;
    private JTable table;
    private JScrollPane scrollPane;
    private LineNumberTable lineTable;
    private JTextField samplename;
    private String idir = "", odir = "";
    protected HelioPlot helioplot;
    protected Data data;
    private HelpSet hs;
    private HelpBroker hb;
    private URL hsURL;
    private ExcelAdapter clipboard;
    protected JRadioButtonMenuItem  uthhebutton, uthsmhebutton,
                                  ternarybutton, logratiobutton;
    private final int N1 = 9, N2 = 10;
    private int numcols  = N1;
    protected boolean pasted = false;
    static final boolean DEBUGGINGMODE = false;
    static final String VERSION = "2.1";
}