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
public class SigmoidComparisonMethods {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        Graph g = ea.matlabToGraph("dblp.mat", "graph");
        // Graph g = ea.matlabToGraph("meshes.mat", "tapir");
        // Graph g = ea.matlabToGraph("meshes.mat", "tapir");
        // Graph g = ea.matlabToGraph("polbooks.mat", "graph");
        //  Graph g = ea.matlabToGraph("netScienceGraph.mat", "graph");
        // Graph g = ea.matlabToGraph("twoMoons.mat", "g6");

        //    Graph g = ea.matlabToGraph("adjnoun.mat", "graph");
        //Graph g = ea.matlabToGraph("airflights.mat", "graph");
        //Graph g = ea.matlabToGraph("epsilonGraphs.mat", "epsilonGraph_02");
        //double[][] coord_iso = new Matrix(ea.readMatlabMatrix("metricLearning.mat", "cTransformed")).getArrayCopy();
        //double[][] coord = ea.readMatlabMatrix("tsneAdj40.mat", "mapped");
        // double[][] coord = ea.readMatlabMatrix("stab.mat", "stab");
        //double[][] coord = ea.readMatlabMatrix("minnesota.mat", "min2d");
        //LINE
        double[][] coordLine = ea.readMatlabMatrix("lineDblp.mat", "o22d");
        double[][] coordGraRep = ea.readMatlabMatrix("dblpGraRep.mat", "cc2d");
        double[][] coordNode2vec = ea.readMatlabMatrix("expNode2vec.mat", "dblp");
        double[][] coordDeepWalk = ea.readMatlabMatrix("expDeepWalk.mat", "dblp");
        // double[][] coord = ea.readMatlabMatrix("res.mat", "res");
        // Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
//        double[][] coordSOM = ea.readMatlabMatrix("resultSOM.mat", "result");
//        GraphCompression gk = new GraphCompression(g, coordSOM);
//        System.out.println("costSOM: ");
//        gk.mdlFunctionSimpleSigmoidComparisonMethods();

        //double[][] coordFR = ea.readMatlabMatrix("Airflights_coordinates.mat", "CoordMatrix2D");
        // double[][] coord = ea.readMatlabMatrix("800_9173.mat", "coord");
        //  double[][] labels = ea.readMatlabMatrix("football.mat", "labels");
        //double[][] labels = ea.readMatlabMatrix("labels_airflights.mat", "labels");
        // double[][] labels = ea.readMatlabMatrix("polbooks.mat", "labels");
//        double[][] labels = ea.readMatlabMatrix("labelsDblp.mat", "labels");
////        // double[][] labels = ea.readMatlabMatrix("reducedLabels.mat", "labels");
////////        //  double[][] labels = ea.readMatlabMatrix("adjnoun.mat", "labels");
//        int[] ids = new int[labels.length];
//        for (int i = 0; i < ids.length; i++) {
//            ids[i] = (int) labels[i][0];
//        }
//        double[][] adjDist = ea.readMatlabMatrix("metricGraph.mat", "AdjDist");
//        double[] ad = new double[adjDist[0].length];
//        for (int i = 0; i < ad.length; i++) {
//            ad[i] = adjDist[0][i];
//        }
//        MetricCorrection mc = new MetricCorrection(ad, 300);
//        mc.run();
//        double[] mm = mc.metric;
//        Visualization v = new Visualization(g);
//        GraphCompression gk = new GraphCompression(g);
//        double cost = gk.mdlFunction(ad);
//        System.out.println(cost);
        // v.displayCoord(coord, new Double(cost).toString());
//        Visualization v = new Visualization(g);
//        System.out.println(g.getVertexCount());
        // double[][] coord_iso = v.getCoordinatesMDS();
//        GraphCompression gk_iso = new GraphCompression(g, coord_iso);
//       double cost_iso = gk_iso.mdlFunction();
//       System.out.println(cost_iso);
        //v.displayCoordNew(coord, " ", ids);
        // v.displayCoordNew(coord, " ");
        // v.displayCoord(coord, new Double(cost_iso).toString());
        //DataUtils du = new DataUtils();
        //du.saveAsMatlab(coord, "iso", "iso10d.mat");
        //ea.writeDoubleToMatlab(coord_iso, "isomap");
        //Check
        DataUtils du = new DataUtils();
        double[][] coordLines = du.scaleLargestAxis(coordLine);
        GraphCompression gk = new GraphCompression(g, coordLines);
        //GraphCompression gk = new GraphCompression(g);
//        double costWithoutEmbedding = gk.codingCostNoEmbedding();
//        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        gk.mdlFunctionSimpleSigmoidComparisonMethods();

        double[][] coordGraReps = du.scaleLargestAxis(coordGraRep);
        gk = new GraphCompression(g, coordGraReps);
        gk.mdlFunctionSimpleSigmoidComparisonMethods();

        double[][] coordNode2vecs = du.scaleLargestAxis(coordNode2vec);
        gk = new GraphCompression(g, coordNode2vecs);
        gk.mdlFunctionSimpleSigmoidComparisonMethods();

        double[][] coordDeepWalks = du.scaleLargestAxis(coordDeepWalk);
        gk = new GraphCompression(g, coordDeepWalks);
        gk.mdlFunctionSimpleSigmoidComparisonMethods();

    }

}
