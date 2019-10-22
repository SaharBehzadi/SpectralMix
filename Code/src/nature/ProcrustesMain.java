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
public class ProcrustesMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        DataUtils du = new DataUtils();
        Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");
       //Graph g = ea.matlabToGraph("meshes.mat", "eppstein");
       Procrustes p = new Procrustes(g, 3);
       p.cluster();
       p.partition();
       p.embedd();
       p.mergeClusters(0, 2);
    }
    
}
