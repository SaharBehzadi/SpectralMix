/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author plantc59cs
 */
public class RunMinnesota {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        // TODO code application logic here
//        String filename = args[0];
//        String varName = args[1];
        IO ea = new IO();
        //Graph g = ea.matlabToGraph(filename, varName);
        //Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");

        // Graph g = ea.matlabToGraph("football.mat", "graph");
        // Graph g = ea.matlabToGraph("airflights.mat", "graph");
        // Graph g = ea.matlabToGraph("polbooks.mat", "graph");
         Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //Graph g = ea.matlabToGraph("epsilonGraphs.mat", "epsilonGraph_02");
        // Graph g = ea.matlabToGraph("sphere3d.mat", "graph");
        //  Graph g = ea.matlabToGraph("twoMoons.mat", "g6");
        //   Graph g = ea.matlabToGraph("adjnoun.mat", "graph");
        // Graph g = ea.matlabToGraph("polblogs.mat", "graph"); //is not connected
        //Graph g = ea.matlabToGraph("can_229.mat", "graph");
        //Graph g = ea.matlabToGraph("sphere3d.mat", "graph");
        //Graph g = ea.matlabToGraph(filename, varName);
        GraphClustering gc = new GraphClustering(g, 0.4);
        gc.run();
        gc.getClusterIds();

        WeightedMajorizationClustered wc = new WeightedMajorizationClustered(g, 2, 1, gc.clusterIds, gc.numClusters);
        wc.run();
        double[][] coord = wc.coord;
        DataUtils du = new DataUtils();
        du.saveAsMatlab(coord, "coord", "result.mat");
//        Visualization v = new Visualization(g);
//        v.displayCoordNew(coord, "bla");

//           
//          for(int i = 0; i < 10; i++){ 
//        HierarchicalEmbeddingWeighted he = new HierarchicalEmbeddingWeighted(g, 0.4);
//            
        //HierarchicalEmbedding he = new HierarchicalEmbedding(g, 0.8);
        // he.run();
        // }
    }

}
