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
public class ClustererMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
        IO ea = new IO();  
        Graph g = ea.matlabToGraph("graph.mat", "graph");
        GraphClusterer gc = new GraphClusterer(g);
        gc.clusterJungEdgeBetweenness(20);
        
    }
}
