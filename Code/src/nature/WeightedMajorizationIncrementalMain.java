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
public class WeightedMajorizationIncrementalMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        DataUtils du = new DataUtils();
        Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");
        WeightedMajorizationIncremental wi = new WeightedMajorizationIncremental(g);
        wi.initialize();
        wi.updateAll();
        double[][] coordIncr = wi.coordNew;
        ea.writeDoubleToMatlab(coordIncr, "incremental");
        
        WeightedMajorizationPlain wp = new WeightedMajorizationPlain();
        double[][] coordPlain = wp.initAndTwoInterations(g);
        ea.writeDoubleToMatlab(coordPlain, "plain");
        
    }

}
