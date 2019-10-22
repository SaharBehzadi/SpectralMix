/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nature;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

/**
 *
 * @author claudia
 */
public class GeneratorMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        //GraphGenerator g = new GraphGenerator(2, 300);
         GraphGenerator gen = new GraphGenerator();
        //Graph<Integer,Integer> gg = g.generateKleinberg();
        //Graph<Integer,Integer> gg = g.generateKleinbergSelf();
        Graph<Integer,Integer> gg = gen.generateThresholdGraph(0.1, 300);
        DataUtils du = new DataUtils();
        Graph g = du.getLargestComponent(gg);
        System.out.println("numVertices: " + g.getVertexCount());
        IO ea = new IO();
        ea.writeGraphToMatlab(g, "g300");
        System.out.println("done");
        

    }

}
