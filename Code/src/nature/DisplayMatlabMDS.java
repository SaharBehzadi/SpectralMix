/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;
import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author claudia
 */
public class DisplayMatlabMDS {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        Graph g = ea.matlabToGraph("metricGraph.mat", "graph");
        
        //double[][] coord_iso = new Matrix(ea.readMatlabMatrix("metricLearning.mat", "cTransformed")).getArrayCopy();
        
           double[][] coord_iso = ea.readMatlabMatrix("coord.mat", "test");


//        double[][] adjDist = ea.readMatlabMatrix("metricGraph.mat", "AdjDist");
//        double[] ad = new double[adjDist[0].length];
//        for (int i = 0; i < ad.length; i++) {
//            ad[i] = adjDist[0][i];
//        }
//        MetricCorrection mc = new MetricCorrection(ad, 300);
//        mc.run();
//        double[] mm = mc.metric;
//        Visualization v = new Visualization(g);
//        GraphCompression gk = new GraphCompression(g);
//        double cost = gk.mdlFunction(ad);
//        System.out.println(cost);


        // v.displayCoord(coord, new Double(cost).toString());

        Visualization v = new Visualization(g);
        
       // double[][] coord_iso = v.getCoordinatesMDS();
        GraphCompression gk_iso = new GraphCompression(g, coord_iso);
       double cost_iso = gk_iso.mdlFunction();
       System.out.println(cost_iso);
       
       v.displayCoord(coord_iso, "iso");
      // v.displayCoord(coord, new Double(cost_iso).toString());
       
        //DataUtils du = new DataUtils();
        //du.saveAsMatlab(coord, "iso", "iso10d.mat");
        //ea.writeDoubleToMatlab(coord_iso, "isomap");
    }
}
