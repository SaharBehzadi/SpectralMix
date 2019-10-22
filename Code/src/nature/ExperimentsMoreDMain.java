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
public class ExperimentsMoreDMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        int dim = 10;
        String filename = args[0];
        String varName = args[1];
        IO ea = new IO();
        Graph g = ea.matlabToGraph(filename, varName);
        
        
        //init with 2D then more D
       // HierarchicalEmbeddingWeighted he = new HierarchicalEmbeddingWeighted(g, 0.4);
        //double[][] coord2D = he.returnCoord();
//        DataUtils du = new DataUtils();
//        double[][] coord2DScaled = du.scaleAndNormalize(coord2D);
//        System.out.println("2D finished");
//        double[][] initMoreD = new double[coord2D.length][dim];
//        Random r = new Random(1);
//        for (int i = 0; i < initMoreD.length; i++) {
//            for (int j = 0; j < initMoreD[i].length; j++) {
//                if (j < 2) {
//                    initMoreD[i][j] = coord2DScaled[i][j];
//                } else {
//                    initMoreD[i][j] = r.nextDouble();
//                }
//            }
//        }
//        WeightedMajorization wj = new WeightedMajorization(g, dim, 20, initMoreD);
//        wj.run(1e-3);
//        double[][] res = wj.bestCoord;
        //init with 2D
        
        //directly moreD
      HierarchicalEmbeddingWeightedMoreD he = new HierarchicalEmbeddingWeightedMoreD(g, 0.4, dim);
      double[][] res = he.run();
      
//      
      DataUtils du = new DataUtils();
       du.saveAsMatlab(res, "coord", "resultMoreD.mat");
//        System.out.println("MoreD finished");
//        double[][] labels = ea.readMatlabMatrix("football.mat", "labels");
//        int[] ids = new int[labels.length];
//        for (int i = 0; i < ids.length; i++) {
//            ids[i] = (int) labels[i][0];
//        }
//        du.writeArff(res, ids, 12, "football");
//      System.out.println("done");
        
    }
}
