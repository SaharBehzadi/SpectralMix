/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;
import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author claudia.plant
 */
public class SimAnnMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        Graph g = ea.matlabToGraph("metricGraph.mat", "graph");
        double[][] coord_groundTruth = ea.readMatlabMatrix("coord.mat", "test");
        
        Visualization v = new Visualization(g);
        double[][] coord_mj = v.getCoordinatesItMajWithInit(coord_groundTruth);
         //double[][] coord_mj = v.getCoordinatesItMaj(coord_groundTruth);
    
//        GraphCompression gk_iso = new GraphCompression(g, coord_mj);
//        double cost_iso = gk_iso.mdlFunction();
//        double[] costs = gk_iso.edgeCosts;
//        System.out.println(cost_iso);

       // v.displayCoord(coord_mj, gk_iso.edgeCosts, "mj");
    }
}
