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
public class GridMain {

    static double convConstant = 1E-3;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        //Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //   Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");
        Graph g = ea.matlabToGraph("dti.mat", "graph"); //random 70 unverfaltet; ist es durch Grid schlechter? 
        double[][] coord = ea.readMatlabMatrix("init.mat", "init");
        double maxS = 20.0;
//        ////////////2D////////////////////////////////////////////////////////////////
        Grid gg = new Grid(g, maxS);
        //gg.run(convConstant);
        gg.run(coord, convConstant);

        double[][] cRes = gg.getCoord();

        ///////////////////2D////////////////////////////////////////////////////////
        ////////////////3D//////////////////////////////////////////////////////////
//        DisplayEmbedding3d dd = new DisplayEmbedding3d(g, coord);
//        dd.display();
//        double[][] cRes = dd.getCoord();
        //3D///////////////////////////////////////////////////////////////////////////////
        Visualization v = new Visualization(g);        
        v.displayCoordNew(cRes, " ");
//        DataUtils du = new DataUtils();
//        du.saveAsMatlab(cRes, "coord", "result.mat");

    }

}
