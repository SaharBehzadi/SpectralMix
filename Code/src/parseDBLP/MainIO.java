/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package parseDBLP;

import edu.uci.ics.jung.graph.Graph;
import nature.IO;

/**
 *
 * @author claudia
 */
public class MainIO {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String[] names = new String[3];
        names[0] = "Christos Faloutsos";
        names[1] = "Philip S. Yu";
        names[2] = "Chris H. Q. Ding";
        Parser pp = new Parser("dblp.xml");
        Person p = new Person();
        Graph g = p.getCoauthors(names);
        int[] ids = p.clid;
        String[] nn = p.nn;
        

        System.out.println(g.getVertexCount() + " vertices " + g.getEdgeCount() + " edges");
        IO ea = new IO();
        ea.writeGraphToMatlab(g, "adj");
    }
}
