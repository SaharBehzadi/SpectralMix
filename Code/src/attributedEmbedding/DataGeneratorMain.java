/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package attributedEmbedding;

import edu.uci.ics.jung.graph.Graph;
import nature.Visualization;

/**
 *
 * @author plantc59cs
 */
public class DataGeneratorMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        DataGenerator dg = new DataGenerator();
        dg.agreementBetweenGraphAndAttributes();
        Graph g = dg.g;
//        Visualization v = new Visualization(g);
//        double[][] coord = v.getCoordinatesIsomapOnly();
//        // double[][] coord = v.getCoordinatesFR();
//        v.displayCoordNew(coord, "classID", dg.classID);
//        v.displayCoordNew(coord, "attribute", dg.attribute);




    }
    
}
