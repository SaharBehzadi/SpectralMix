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
public class InitializationPhaseMain {

    /**
     * @param args the command line arguments
     */
     public static void main(String[] args) throws InterruptedException {
        // TODO code application logic here
        IO ea = new IO();
        Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //Graph g = ea.matlabToGraph("meshes.mat", "smallmesh"); //random 70 unverfaltet; ist es durch Grid schlechter? 
        System.out.println("numEdges: " + g.getEdgeCount());
        int numT = 6;
        int d = 2;
        Parallel gp = new Parallel(g, numT, d);
//        //double[][] init = gp.getBestInitialization(10, 1000, 3);
//        DataUtils du = new DataUtils();
//        du.saveAsMatlab(init, "init", "init.mat");
        
    }
}
