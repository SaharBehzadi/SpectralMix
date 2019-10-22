/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.algorithms.cluster.EdgeBetweennessClusterer;
import edu.uci.ics.jung.graph.Graph;
import java.util.Iterator;
import java.util.Set;

/**
 *
 * @author claudia.plant
 */
public class GraphClusterer {

    Graph g;
    double[] clid;

    public GraphClusterer(Graph g) {
        this.g = g;
        clid = new double[g.getVertexCount()];
    }

    public int clusterJungEdgeBetweenness(int numEdgesToRemove) {
        EdgeBetweennessClusterer eb = new EdgeBetweennessClusterer(numEdgesToRemove);
        Set<Set<Integer>> res = eb.transform(g);
        int clCount = 0;
        for (Iterator<Set<Integer>> cIt = res.iterator(); cIt.hasNext();) {
            Set<Integer> vertices = cIt.next();
            Object[] vv = vertices.toArray();
            for (int i = 0; i < vv.length; i++) {
                int obj = ((Integer) vv[i]).intValue();
                clid[obj] = clCount;

            }
            clCount++;
        }
       return clCount;

    }
}
