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
public class EvaluationRepresentationLearningMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        DataUtils du = new DataUtils();
        IO ea = new IO();
        Graph g = du.readAdjList("dblpAdj");
        double[][] labels = ea.readMatlabMatrix("labelsDblp.mat", "labels");
        //double[][] labels = ea.readMatlabMatrix("labels3.mat", "labels3");
        int[] id = new int[labels.length];
        for (int i = 0; i < id.length; i++) {
            id[i] = (int) labels[i][0];
        }
//        //Gempe
        double[][] coordGempe = du.readCoordFile("coord", g.getVertexCount(), 2);
        EvaluationNodeClassification ev = new EvaluationNodeClassification(g, coordGempe, id, 9, 2);
        double[][] resGempe = ev.nodeClassificationExperiments();
        du.saveAsMatlab(resGempe, "resGempe", "resGempe.mat");
        
        //Verse
        double[][] coordVerse = du.readCoordFile("dblpVerseTsne", g.getVertexCount(), 2);
        ev = new EvaluationNodeClassification(g, coordVerse, id, 9, 2);
        double[][] resVerse = ev.nodeClassificationExperiments();
        du.saveAsMatlab(resVerse, "resVerse", "resVerse.mat");
        
        //node2vec
        double[][] coordNode2vec = du.readCoordFile("dblpNode2vecTsne", g.getVertexCount(), 2);
        ev = new EvaluationNodeClassification(g, coordNode2vec, id, 9, 2);
        double[][] resNode2vec = ev.nodeClassificationExperiments();
        du.saveAsMatlab(resNode2vec, "resNode2vec", "resNode2vec.mat");
        
        //deepWalk
        double[][] coordDeepWalk = du.readCoordFile("dblpDeepWalkTsne", g.getVertexCount(), 2);
        ev = new EvaluationNodeClassification(g, coordDeepWalk, id, 9, 2);
        double[][] resDeepWalk = ev.nodeClassificationExperiments();
        du.saveAsMatlab(resDeepWalk, "resDeepWalk", "resDeepWalk.mat");
        
        //line
        double[][] coordLine = du.readCoordFile("dblpLineTsne", g.getVertexCount(), 2);
        ev = new EvaluationNodeClassification(g, coordLine, id, 9, 2);
        double[][] resLine = ev.nodeClassificationExperiments();
        du.saveAsMatlab(resLine, "resLine", "resLine.mat");
        
//GraRep
 //line
        double[][] coordEnt = du.readCoordFile("dblpAdj-maxent.coord", g.getVertexCount(), 2);
        ev = new EvaluationNodeClassification(g, coordEnt, id, 9, 2);
        double[][] resEnt = ev.nodeClassificationExperiments();
        du.saveAsMatlab(resEnt, "resEnt", "resEnt.mat");
        

    }

}
