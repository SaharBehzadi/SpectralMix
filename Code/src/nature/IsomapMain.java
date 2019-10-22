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
public class IsomapMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
     // Graph g = ea.matlabToGraph("airflights.mat", "graph");
       // Graph g = ea.matlabToGraph("football.mat", "graph");
       //  Graph g = ea.matlabToGraph("sphere3d.mat", "graph");
         //         Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
      //  Graph g = ea.matlabToGraph("polbooks.mat", "graph");
        
        //Graph g = ea.matlabToGraph("graph.mat", "epsilonGraph_02");
            Graph g = ea.matlabToGraph("meshes.mat", "tapir");
               //Graph g = ea.matlabToGraph("graph.mat", "graph");
           //Graph g = ea.matlabToGraph("twoMoons.mat", "g6");
               // Graph g = ea.matlabToGraph("can_229.mat", "graph");
       
        Visualization v = new Visualization(g);
        double[][] coord = v.getCoordinatesIsomapOnly();
       // double[][] coord = v.getLaplacianCoord(2);
        
        //double[][] labels = ea.readMatlabMatrix("football.mat", "labels");
        //double[][] labels = ea.readMatlabMatrix("polbooks.mat", "labels");
//         double[][] labels = ea.readMatlabMatrix("labels_airflights.mat", "labels");
//        int[] ids = new int[labels.length];
//        for (int i = 0; i < ids.length; i++) {
//            ids[i] = (int) labels[i][0];
//        }

       // v.displayCoordNew(coord, " ", ids);
          v.displayCoordNew(coord, " ");
         //   v.displayDegreeDistribution(coord, " ");
        
        DataUtils du = new DataUtils();
        double[][] coords = du.scaleLargestAxis(coord);
        //du.saveAsMatlab(coord, "result", "laplacian.mat");
        GraphCompression gk = new GraphCompression(g, coords);
         double costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
         gk.mdlFunctionSimpleSigmoidComparisonMethods();
        //System.out.println("compression cost: " + compressionCost);

    }
}
