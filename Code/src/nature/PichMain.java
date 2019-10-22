/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;
import edu.uci.ics.jung.graph.Graph;
import java.util.Random;
import mdsj.MDSJ;

/**
 *
 * @author claudia.plant
 */
public class PichMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
       // Graph g = ea.matlabToGraph("airflights.mat", "graph");
        //Graph g = ea.matlabToGraph("football.mat", "graph");
          //Graph g = ea.matlabToGraph("sphere3d.mat", "graph");
           //       Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //Graph g = ea.matlabToGraph("polbooks.mat", "graph");
        
         //Graph g = ea.matlabToGraph("epsilonGraphs.mat", "epsilonGraph_02");
           Graph g = ea.matlabToGraph("twoMoons.mat", "g6");
               // Graph g = ea.matlabToGraph("can_229.mat", "graph");
       
        Visualization v = new Visualization(g);
        //double[][] coord = v.getCoordinatesPich();
        double[][] init = new double[g.getVertexCount()][2];
        Random r = new Random(1);
        for(int i = 0; i < init.length; i++)
            for(int j = 0; j < init[i].length; j++)
                init[i][j] = r.nextDouble();
        
        
        int numObj = g.getVertexCount();
         double[][] dist = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j) || (i == j)) {
                    dist[i][j] = 0.0;
                } else {
                    dist[i][j] = 1.0;
                }
            }
        }
        init = MDSJ.classicalScaling(dist); // apply MDS
        init = new Matrix(init).transpose().getArrayCopy();
        
        double[][] coord = v.getCoordinatesStress(init);
        
//        double[][] labels = ea.readMatlabMatrix("football.mat", "labels");
//        //double[][] labels = ea.readMatlabMatrix("polbooks.mat", "labels");
//         //  double[][] labels = ea.readMatlabMatrix("labels_airflights.mat", "labels");
//        int[] ids = new int[labels.length];
//        for (int i = 0; i < ids.length; i++) {
//            ids[i] = (int) labels[i][0];
//        }

       // v.displayCoordNew(coord, " ", ids);
//       Matrix cm = new Matrix(coord).transpose();
//       coord = cm.getArrayCopy();
          v.displayCoordNew(coord, " ");
        
        DataUtils du = new DataUtils();
        double[][] coords = du.scaleLargestAxis(coord);
        du.saveAsMatlab(coords, "result", "resultNonMetric.mat");
        GraphCompression gk = new GraphCompression(g, coords);
         double costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
         gk.mdlFunctionSimpleSigmoidComparisonMethods();
        //System.out.println("compression cost: " + compressionCost);

    }
}
