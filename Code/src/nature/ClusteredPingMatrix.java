/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import java.util.Random;

/**
 *
 * @author claudia.plant: generates .arrf File with the corresponding
 * categorical variables
 */
public class ClusteredPingMatrix {
    Graph g;
    int numPings;
    boolean verbose = true;

    public ClusteredPingMatrix(Graph g, int numPings) {
        this.g = g;
        this.numPings = numPings;
    }

    public void writePingMatrixToMatlab() {
        double[][] pingMatrix = new double[numPings][g.getVertexCount()];
        Random r = new Random(1);
           System.out.println(g.getEdgeCount());
        for (int i = 0; i < pingMatrix.length; i++) {
             GraphClusterer gc = new GraphClusterer(g);
            int nCl = gc.clusterJungEdgeBetweenness(100 + r.nextInt(20));
        
            while (nCl == 1) {
                gc = new GraphClusterer(g);
                nCl = gc.clusterJungEdgeBetweenness(100 + r.nextInt(20));
            }
            if(verbose)
                System.out.println("Ping " +  i + " : " + nCl + " clusters.");
            pingMatrix[i] = gc.clid;
        }
        IO ea = new IO();
        ea.writeDoubleToMatlab(pingMatrix, "pingMatrix", "pingMatrix");
    }
}
