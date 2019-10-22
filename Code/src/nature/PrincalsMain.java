/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import weka.core.Instances;

/**
 *
 * @author claudia.plant
 */
public class PrincalsMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        //consider each edge as an own variable
           Graph g = ea.matlabToGraph("graph.mat", "graph");
//            Instances dd = ea.graphToArff(g);
//            Princals p = new Princals(dd, 2);
//            p.run();
           int[] ids = new int[g.getVertexCount()];
          
        
        //19.07. Use clusterings from Jung as Pings
        ArffFileReader af = new ArffFileReader();
        Instances dd = af.readFile("clPings.arff");
        Princals p = new Princals(dd, 2);
        p.run();
        double[][] coord = p.objectCoord.transpose().getArrayCopy();
        Visualization vv = new Visualization(g);
        vv.displayCoordSmall(coord, "princals", ids);
        double[][] coordt = p.objectCoord.getArrayCopy();
        GraphCompression gk = new GraphCompression(g, coordt);
        System.out.println(gk.mdlFunction());
        
    }
}
