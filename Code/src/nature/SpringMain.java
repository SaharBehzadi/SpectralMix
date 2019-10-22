package nature;

import edu.uci.ics.jung.graph.Graph;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author claudia.plant
 */
public class SpringMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
         // TODO code application logic here
         IO ea = new IO();
       // Graph g = ea.matlabToGraph("airflights.mat", "graph");
        //Graph g = ea.matlabToGraph("football.mat", "graph");
         Graph g = ea.matlabToGraph("polbooks.mat", "graph");
        Visualization v = new Visualization(g);
        double[][] coord = v.getCoordinatesSpring();
         //double[][] labels = ea.readMatlabMatrix("football.mat", "labels");
          double[][] labels = ea.readMatlabMatrix("polbooks.mat", "labels");
        int[] ids = new int[labels.length];
        for (int i = 0; i < ids.length; i++) {
            ids[i] = (int) labels[i][0];
        }
                v.displayCoordNew(coord, " ", ids);
    
    }
}
