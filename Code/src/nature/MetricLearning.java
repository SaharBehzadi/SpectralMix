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
public class MetricLearning {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        DataUtils du = new DataUtils();
//        double[][] dist2d = du.readMatlabMatrix("metricLearning.mat", "iso2dDist");
//        double[] d2d = new double[dist2d[0].length];
//        for (int i = 0; i < dist2d[0].length; i++) {
//            d2d[i] = dist2d[0][i];
//        }
//
//
//        double[][] dist10d = du.readMatlabMatrix("metricLearning.mat", "iso10dDist");
//        double[] d10d = new double[dist10d[0].length];
//        for (int i = 0; i < dist10d[0].length; i++) {
//            d10d[i] = dist10d[0][i];
//        }
//
//        double[][] dist20d = du.readMatlabMatrix("metricLearning.mat", "iso20dDist");
//        double[] d20d = new double[dist20d[0].length];
//        for (int i = 0; i < dist20d[0].length; i++) {
//            d20d[i] = dist20d[0][i];
//        }
//
//        double[][] dist100d = du.readMatlabMatrix("metricLearning.mat", "iso100dDist");
//        double[] d100d = new double[dist100d[0].length];
//        for (int i = 0; i < dist100d[0].length; i++) {
//            d100d[i] = dist100d[0][i];
//        }
//
//        double[][] dist160d = du.readMatlabMatrix("metricLearning.mat", "iso160dDist");
//        double[] d160d = new double[dist160d[0].length];
//        for (int i = 0; i < dist160d[0].length; i++) {
//            d160d[i] = dist160d[0][i];
//        }

        double[][] pDist = du.readMatlabMatrix("metricLearning.mat", "distRotated");
        double[] dp = new double[pDist[0].length];
        for (int i = 0; i < pDist[0].length; i++) {
            dp[i] = pDist[0][i];
        }


        IO ea = new IO();
        Graph g = ea.matlabToGraph("metricGraph.mat", "graph");

        GraphCompression gk = new GraphCompression(g);
        //double cost2d = gk.mdlFunction(d2d);
        //gk.writeContraintsMatrix(d10d);
        // double cost10d = gk.mdlFunction(d10d);
//        double cost20d = gk.mdlFunction(d20d);
//        double cost100d = gk.mdlFunction(d100d);
//        double cost160d = gk.mdlFunction(d160d);
        double costPath = gk.mdlFunction(dp);
        gk.writeContraintsMatrix(dp);
        // System.out.println("2d: " + cost2d);
        //System.out.println("10d: " + cost10d);
//        System.out.println("20d: " + cost20d);
//        System.out.println("100d: " + cost100d);
//        System.out.println("160d: " + cost160d);
        System.out.println("rotated " + costPath);



    }
}
