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
public class GraphClusteringMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");
        GraphClustering gc = new GraphClustering(g, 0.3);
        gc.run();
        gc.getClusterIds();
        int[][] ids = gc.clusterIds;
        System.out.println("m");

    }
}
