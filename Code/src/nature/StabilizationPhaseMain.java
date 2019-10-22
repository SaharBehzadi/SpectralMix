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
public class StabilizationPhaseMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //Graph g = ea.matlabToGraph("meshes.mat", "smallmesh"); //random 70 unverfaltet; ist es durch Grid schlechter? 
        System.out.println("numEdges: " + g.getEdgeCount());
        DataUtils du = new DataUtils();
        double[][] coordInit = du.readMatlabMatrix("init.mat", "init");
        int numT = 6;
        Parallel gp = new Parallel(g, numT, coordInit, true);
//        double[][] res = gp.stabilizationPhase(1000, 3);
//        du.saveAsMatlab(res, "stabilization", "stab.mat");
    }
    
}
