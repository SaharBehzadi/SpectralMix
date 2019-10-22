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
public class EvaluationNodeClusteringMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        DataUtils du = new DataUtils();
         Graph g = du.readAdjList("pubmedAdj");
        double[][] labels = ea.readMatlabMatrix("labelsDiabetes.mat", "labels");
        //double[][] labels = ea.readMatlabMatrix("labels3.mat", "labels3");
        int d = 128;
        int numClasses = 3;
        int[] id = new int[labels.length];
        for (int i = 0; i < id.length; i++) {
            id[i] = (int) labels[i][0];
        }
//        //Gempe
        double[][] coordGempe = du.readCoordFile("pubmedLine.txt", g.getVertexCount(), d);
        EvaluationNodeClustering ev = new EvaluationNodeClustering(g, coordGempe, id, numClasses, d);
        ev.nodeClustering();
    }
    
}
