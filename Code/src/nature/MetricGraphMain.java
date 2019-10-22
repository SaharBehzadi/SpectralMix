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
public class MetricGraphMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        
        //Graph g = ea.matlabToGraph("metricGraph.mat", "graph");
        
        DataUtils du = new DataUtils();
        double[][] dd = du.readMatlabMatrix("metricGraph.mat", "dist");
        double[] dist = new double[dd[0].length];
        for (int i = 0; i < dd[0].length; i++){
            if(dd[0][i] == 0) //no edge
                dist[i] = 100;
            if(dd[0][i] == 1) //edge
                dist[i] = 1;
            //dist[i] = dd[0][i];
        }
        MetricCorrection mc = new MetricCorrection(dist, 300);
        mc.run();
        double[] mm = mc.metric;
        double[][] res = new double[mm.length][1];
        for(int i = 0; i < mm.length; i++)
            res[i][0] = mm[i];
        
        du.saveAsMatlab(res, "res", "distance.mat");
        
        
//        Visualization v = new Visualization(g);
//         //double[][] coord = v.getCoordinatesClusteredPCA(1000, 2);
//        double[][] coord = v.getCoordinatesIsomapOnly();
//        coord = new Matrix(coord).transpose().getArrayCopy();
//      
//        
//          // double[][] coord = v.getCoordinatesPingICANeighbors(2, 1000);
//           
//        //parameters numSteps, numWalks
////        double[][] coord = v.getCoordinatesPingPCA(50, 2);
//        GraphCompression gk = new GraphCompression(g, coord);
//        System.out.println(gk.mdlFunction());
//        v.displayCoord(coord, "Isomap");
        
    }
}
