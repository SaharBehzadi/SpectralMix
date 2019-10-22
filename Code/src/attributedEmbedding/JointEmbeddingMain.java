/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package attributedEmbedding;

import datautils.Sociopatterns;
import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author plantc59cs
 */
public class JointEmbeddingMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        // TODO code application logic here
        DataGenerator dg = new DataGenerator();
       dg.agreementBetweenGraphAndAttributes();
//       //dg.randomGraphRandomAttributes();
//        //dg.disagreementBetweenGraphAndAttributes();
//        //dg.clusteredGraphRandomAttributes();
//        //dg.randomGraphRandomAttributes();
Graph g = dg.g;
//        Visualization v = new Visualization(g);
//        double[][] coord = v.getCoordinatesIsomapOnly();
//        // double[][] coord = v.getCoordinatesFR();
//        v.displayCoordNew(coord, "classID", dg.classID);
//        v.displayCoordNew(coord, "attribute", dg.attribute);

//        Sociopatterns sp = new Sociopatterns();
//        sp.readData();
//        Graph g = sp.getMeasured();

//
        int[] classId = new int[g.getVertexCount()];
        for (int i = 0; i < classId.length; i++) {
            if (i > 500) {
                classId[i] = 1;
            }
        }
        int numEdges = g.getEdgeCount();
        int numCat = g.getVertexCount();
  double factor = (double) numEdges / (double) numCat;
//        int numEdgesGraph = sp.getContactsSumOfWeights();
//        int numEdgesCat = sp.getGenderSumOfWeights();
//        double factor = (double) numEdgesGraph / (double) numEdgesCat;
//check
  //      double dd = numEdgesCat * factor;
        //JointEmbedding je = new JointEmbedding(g, sp.getGender(), 2, 3, sp.getClassId(), 1.0, factor);
         JointEmbedding je = new JointEmbedding(g, dg.attribute, 2, 3, dg.classID, 1.0, factor);
        je.init(0);
        je.runWeighted();
//        SpectralEmbedding se = new SpectralEmbedding(g, 3, dg.classID);
//        se.init(0);
//        se.run();
//        se.checkEigenvector();
//        System.out.println("m");
    }

}
