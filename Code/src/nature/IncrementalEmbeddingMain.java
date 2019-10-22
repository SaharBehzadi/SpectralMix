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
public class IncrementalEmbeddingMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
         IO ea = new IO();
        //Graph g = ea.matlabToGraph(filename, varName);
       //Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");
         //Graph g = ea.matlabToGraph("polbooks.mat", "graph");
         // Graph g = ea.matlabToGraph("football.mat", "graph");
      //Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        Graph g = ea.matlabToGraph("graph.mat", "graph");
      //Graph g = ea.matlabToGraph("airflights.mat", "graph");
       IncrementalEmbedding ie = new IncrementalEmbedding(g);
       ie.run();
    }
    
}
