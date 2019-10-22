/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author claudia.plant
 */
public class PingMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
//         IO ea = new IO();
//        double[][] adj = ea.readMatlabMatrix("graph.mat", "graph");
//        Ping p = new Ping(adj);
//        double[][] pMatrix = p.pingMatrix();
//        ea.writeDoubleToMatlab(pMatrix, "pMatrix.mat");
//        
        
        
        IO ea = new IO();
        
        Graph g = ea.matlabToGraph("graph.mat", "graph");
        Visualization v = new Visualization(g);
         //double[][] coord = v.getCoordinatesClusteredPCA(1000, 2);
           double[][] coord = v.getCoordinatesPingICANeighbors(2, 1000);
           
        //parameters numSteps, numWalks
//        double[][] coord = v.getCoordinatesPingPCA(50, 2);
        GraphCompression gk = new GraphCompression(g, coord);
        System.out.println(gk.mdlFunction());
      
        
    }
}
