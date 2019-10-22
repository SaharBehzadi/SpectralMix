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
public class IncrementalParallelMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        DataUtils du = new DataUtils();
        //Graph g = du.readAdjList("dblpAdj");
        //double[][] init = du.readMatlabMatrix("initMinn100.mat", "init");
        //double[][] init = du.readMatlabMatrix("init.mat", "init");
        //Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //Graph g = ea.matlabToGraph("football.mat", "graph");
        //Graph g = ea.matlabToGraph("meshes.mat", "eppstein");
        //Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");
        //Graph g = ea.matlabToGraph("bus.mat", "graph");
        Graph g = ea.matlabToGraph("powerGrid.mat", "graph");
        //Graph g = du.readLuxemburg();
        //double[][] coord = ea.readMatlabMatrix("lCoord.mat", "coord");
        //double[][] coord = ea.readMatlabMatrix("initMinn100.mat", "init");
        //double[][] coord = ea.readMatlabMatrix("init.mat", "init");
        //double[][] labels = ea.readMatlabMatrix("football.mat", "labels");
//         int[] id = new int[labels.length];
//         for(int i = 0; i < id.length; i++)
//             id[i] = (int)labels[i][0];
//         Visualization v = new Visualization(g);
//        // v.displayCoordNew(coord, "bla", id);
//         v.displayCoordNew(coord, "bla");

        //System.out.println("m");
        IncrementalParallel p = new IncrementalParallel(g);

        //public void getBestInitalization(int numTry, int numUpdatesPerEdge, int numT)
        p.getBestInitalization(10, 10, 8);
        //public void refine(int numUpdatesPerEdge, int numT, int numClusters) {
        //p.refine(1, 1, 3);

        IncrementalParallelGrid pg = new IncrementalParallelGrid(g, p.getBestCoord());
        //IncrementalParallelGrid pg = new IncrementalParallelGrid(g, init);
        // public void run(int iter, int numUpdatesPerEdge,  int numT)
        pg.run(10, 1, 8);
        //p.run(2, 1000);
    }

}
