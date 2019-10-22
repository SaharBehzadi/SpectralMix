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
public class AntrittsvorlesungMain {

    public static void main(String[] args) {
        IO ea = new IO();
        //Graph g = ea.matlabToGraph(filename, varName);
        //Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");
         Graph g = ea.matlabToGraph("similarityGraph.mat", "similarity");
         System.out.println(g.getVertexCount() + " " + g.getEdgeCount());
//        GraphGenerator gen = new GraphGenerator();
//        Graph g = gen.generateEmptyGraph(200);

         double[][] data = ea.readMatlabMatrix("twoMoons500.mat", "coord");
//        //Graph g = ea.matlabToGraph("football.mat", "graph");
//        double[][] labels = ea.readMatlabMatrix("dataFigNew.mat", "labelsNat");
//        //double[][] labels = ea.readMatlabMatrix("dataFigNew.mat", "dummy");
//        double[][] coordCl = ea.readMatlabMatrix("dataFigNew.mat", "coordNew");
//        double[][] coordU = ea.readMatlabMatrix("dataFigNew.mat", "cu");
        //double[][] coord = ea.readMatlabMatrix("dataFigNew.mat", "rr");

        //double[][] labels = ea.readMatlabMatrix("labels_airflights.mat", "labels");
        //double[][] labels = ea.readMatlabMatrix("polbooks.mat", "labels");
//        double[][] labels = ea.readMatlabMatrix("labelsDblp.mat", "labels");
//        // double[][] labels = ea.readMatlabMatrix("reducedLabels.mat", "labels");
//////        //  double[][] labels = ea.readMatlabMatrix("adjnoun.mat", "labels");
//        int[] ids = new int[labels.length];
//        for (int i = 0; i < ids.length; i++) {
//            ids[i] = (int) labels[i][0];
//        }
      //  Graph g = gen.generateClusteredGraph(ids, coordCl);
        Visualization v = new Visualization(g);
        // double[][] b = v.getCoordinatesCircle();
        // v.displayCoordNew(coord, "circle");
        v.displayCoordNew(data, " ");
    }

}
