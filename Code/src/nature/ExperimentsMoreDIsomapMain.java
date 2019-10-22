/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import java.util.Random;

/**
 *
 * @author plantc59cs
 */
public class ExperimentsMoreDIsomapMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        int dim = 5;
        String filename = "airflights.mat";
        String varName = "graph";
        IO ea = new IO();
        Graph g = ea.matlabToGraph(filename, varName);
        Visualization v = new Visualization(g);
//        double[][] coord2D = v.getCoordinatesIsomapOnly();
        DataUtils du = new DataUtils();
//        du.saveAsMatlab(coord2D, "result", "resultIsomap.mat");
//        System.out.println("2D finished");
       
        double[][] labels = ea.readMatlabMatrix("labels_airflights.mat", "labels");
        int[] ids = new int[labels.length];
        for (int i = 0; i < ids.length; i++) {
            ids[i] = (int) labels[i][0];
        }
        
        double[][] coord2D = ea.readMatlabMatrix("resultISOMAP.mat", "result");
        
        du.writeArff(coord2D, ids, 7, "football");
        System.out.println("done");
        
    }
}
