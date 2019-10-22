


package nature;
/*
 * Visu2D.java
 *
 * Created on 07 February 2005, 15:48
 */

/**
 *
 * @author  Administrator
 */

import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import javax.swing.*;



public class VisuInfo extends javax.swing.JFrame {
    DataObject[] clustered;
    int[][] internalCoord;
    String title;
    int dim0;
    int dim1;
    //Ver�ndert f�r Darstellung Beispiel
    Color [] cols = {Color.GREEN, Color.BLUE, Color.RED, Color.LIGHT_GRAY, Color.CYAN, Color.YELLOW, Color.PINK, Color.DARK_GRAY, Color.MAGENTA, Color.getHSBColor(0.0f, 0.5f, 1.0f), Color.getHSBColor(0.1f, 0.5f, 1.0f), Color.getHSBColor(0.2f, 0.5f, 1.0f), Color.getHSBColor(0.3f, 0.5f, 1.0f), Color.getHSBColor(0.4f, 0.5f, 1.0f), Color.getHSBColor(0.5f, 0.5f, 1.0f), Color.getHSBColor(0.6f, 0.5f, 1.0f), Color.getHSBColor(0.7f, 0.5f, 1.0f), Color.getHSBColor(0.8f, 0.5f, 1.0f), Color.getHSBColor(0.9f, 0.5f, 1.0f), Color.getHSBColor(1.0f, 0.5f, 1.0f), Color.GREEN, Color.BLUE, Color.RED, Color.ORANGE, Color.LIGHT_GRAY, Color.CYAN, Color.YELLOW, Color.PINK, Color.DARK_GRAY, Color.MAGENTA, Color.getHSBColor(0.0f, 0.5f, 1.0f), Color.getHSBColor(0.1f, 0.5f, 1.0f), Color.getHSBColor(0.2f, 0.5f, 1.0f), Color.getHSBColor(0.3f, 0.5f, 1.0f), Color.getHSBColor(0.4f, 0.5f, 1.0f), Color.getHSBColor(0.5f, 0.5f, 1.0f), Color.getHSBColor(0.6f, 0.5f, 1.0f), Color.getHSBColor(0.7f, 0.5f, 1.0f), Color.getHSBColor(0.8f, 0.5f, 1.0f), Color.getHSBColor(0.9f, 0.5f, 1.0f), Color.getHSBColor(1.0f, 0.5f, 1.0f)} ;
    //Color [] cols = {Color.GREEN, Color.BLUE, Color.LIGHT_GRAY, Color.ORANGE, Color.RED, Color.CYAN, Color.YELLOW, Color.PINK, Color.DARK_GRAY, Color.MAGENTA, Color.getHSBColor(0.0f, 0.5f, 1.0f), Color.getHSBColor(0.1f, 0.5f, 1.0f), Color.getHSBColor(0.2f, 0.5f, 1.0f), Color.getHSBColor(0.3f, 0.5f, 1.0f), Color.getHSBColor(0.4f, 0.5f, 1.0f), Color.getHSBColor(0.5f, 0.5f, 1.0f), Color.getHSBColor(0.6f, 0.5f, 1.0f), Color.getHSBColor(0.7f, 0.5f, 1.0f), Color.getHSBColor(0.8f, 0.5f, 1.0f), Color.getHSBColor(0.9f, 0.5f, 1.0f), Color.getHSBColor(1.0f, 0.5f, 1.0f), Color.GREEN, Color.BLUE, Color.RED, Color.ORANGE, Color.LIGHT_GRAY, Color.CYAN, Color.YELLOW, Color.PINK, Color.DARK_GRAY, Color.MAGENTA, Color.getHSBColor(0.0f, 0.5f, 1.0f), Color.getHSBColor(0.1f, 0.5f, 1.0f), Color.getHSBColor(0.2f, 0.5f, 1.0f), Color.getHSBColor(0.3f, 0.5f, 1.0f), Color.getHSBColor(0.4f, 0.5f, 1.0f), Color.getHSBColor(0.5f, 0.5f, 1.0f), Color.getHSBColor(0.6f, 0.5f, 1.0f), Color.getHSBColor(0.7f, 0.5f, 1.0f), Color.getHSBColor(0.8f, 0.5f, 1.0f), Color.getHSBColor(0.9f, 0.5f, 1.0f), Color.getHSBColor(1.0f, 0.5f, 1.0f)} ;
    
    // innere Klasse f�r Zeichenfl�che
    private class Leinwand extends JPanel implements MouseListener{
        
        public void mousePressed(MouseEvent e){
            
        }
        public void mouseReleased(MouseEvent e){
            
        }
        
        public void mouseEntered(MouseEvent e){
            
        }
        
        public void mouseExited(MouseEvent e){
            
        }
        
        public void mouseClicked(MouseEvent e){
            int x = e.getX();
            int y = e.getY();
            double minDist = Double.MAX_VALUE;
            int minIndex = clustered.length;
            for(int i = 0; i < internalCoord.length; i++){
                double aktDist = distance(x, y, internalCoord[i][0], internalCoord[i][1]);
                if(aktDist < minDist){
                    minDist = aktDist;
                    minIndex = i;
                }
            }
            //count how many objects are in this location
            int counter = 0;
             for(int i = 0; i < internalCoord.length; i++){
              if (distance(x, y, internalCoord[i][0], internalCoord[i][1]) == minDist){
               counter++;
            }
             }
            String s = clustered[minIndex].getObjectInfo();
            JFrame info = new JFrame();
            info.setTitle("object info " + "(" + title +") " + Integer.toString(counter));
            info.setLocation(x, y);
            info.setSize(500, 150);
            info.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            JLabel displayLabel = new JLabel();
            info.getContentPane().add(displayLabel);
            displayLabel.setText(s);
            info.setVisible(true);
        }
        
        private double distance(int x1, int y1, int x2, int y2){
            return Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
        }
        
        public void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setStroke(new BasicStroke(6.0f));
            //hier zeichnen, Koordinaten normieren
            double max_x = -Double.MAX_VALUE;
            double max_y = -Double.MAX_VALUE;
            double min_x = Double.MAX_VALUE;
            double min_y = Double.MAX_VALUE;
            for(int i = 0; i < clustered.length; i++){
                if(clustered[i].getCoord()[dim0] > max_x){
                    max_x = clustered[i].getCoord()[dim0];
                }
                if(clustered[i].getCoord()[dim1] > max_y){
                    max_y = clustered[i].getCoord()[dim1];
                }
                if(clustered[i].getCoord()[dim0] < min_x){
                    min_x = clustered[i].getCoord()[dim0];
                }
                if(clustered[i].getCoord()[dim1] < min_y){
                    min_y = clustered[i].getCoord()[dim1];
                }
            }
            for(int i = 0; i < clustered.length; i++){
                float xCoord = (float) ((clustered[i].getCoord()[dim0]-min_x)/(max_x-min_x) * 500);
                internalCoord[i][0] = (int)xCoord;
                float yCoord = (float) (500 - (clustered[i].getCoord()[dim1]-min_y)/(max_y-min_y) * 500);
                internalCoord[i][1] = (int)yCoord;
                Line2D.Float akt = new Line2D.Float(xCoord, yCoord, xCoord, yCoord);
                if(clustered[i].classID < 40)
                    g2.setColor(cols[clustered[i].classID]);
//                //TEST
                if(clustered[i].classID == 100)
                    g2.setColor(Color.BLACK);
//                //TEST
                g2.setStroke(new BasicStroke((float)(clustered[i].numericClassID * 10)));
                g2.draw(akt);
            }
        }
    }
    
    
    public VisuInfo(DataObject[] store, String title, int dim0, int dim1){
        internalCoord = new int[store.length][2];
//        if(store[0].getCoord().length > 2){
//            PCA p = new PCA(store);
//            DataObject[] reduced = p.pca(2);
//            this.clustered = reduced;
//        } else{
            this.clustered = store;
            //this.clusters = clusters;
        //}
        // Hauptfenster einrichten
        this.title = title;
        this.dim0 = dim0;
        this.dim1 = dim1;
        setTitle(title);
        Leinwand l = new Leinwand();
        l.addMouseListener(l);
        l.setBackground(Color.WHITE);
        getContentPane().add(l, BorderLayout.CENTER);
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    }
}


