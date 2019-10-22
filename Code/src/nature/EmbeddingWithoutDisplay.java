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
public class EmbeddingWithoutDisplay {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//        String filename = args[0];
//        String varName = args[1];
        // TODO code application logic here
        IO ea = new IO();
        //Graph g = ea.matlabToGraph("graph.mat", "graph");
        Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");
        // Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //  Graph g = ea.matlabToGraph("smallTest.mat", "graph");

        //Graph g = ea.matlabToGraph(filename, varName);
        //Graph g = ea.matlabToGraph("meshes.mat", "eppstein");

//        double[][] labels = ea.readMatlabMatrix("labels.mat", "labels");
//        int[] l = new int[labels.length];
//        for(int i  = 0; i < l.length; i++)
//            l[i] = (int) labels[i][0];
//        

//        double[][] coordInit = ea.readMatlabMatrix("resultOwn.mat", "coord");
//
        Visualization v = new Visualization(g);
//
        HierarchicalWeightedGrid gg = new HierarchicalWeightedGrid(g, 20);
        gg.initialize();
        gg.run();
//        
//        WeightedMajorization wj = new WeightedMajorization(g, 2, 10, 1.0);
//        double coord[][] = wj.run();

        double bestCost = gg.getBestCost();
        double savedBits = gg.getSavedBits();
        int bestIter = gg.getBestIteration();
        double[][] coord = gg.getBestDb();

//        double bestCost = wj.bestCost;
//        double savedBits = wj.savedBits;
//        int bestIter = wj.bestIteration;
//        




        System.out.println("Best cost: " + bestCost + " savedBits: " + savedBits + " in iter " + bestIter);
//
//        ea.writeDoubleToMatlab(coord, "result", "result2d");
//        // System.out.println("m");
        v.displayCoordNew(coord, "bla");
    }
}
