/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author claudia.plant
 */
public class IOMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        //ea.pingMatrixToArff("pingMatrix.mat", "pingMatrix");
        ea.writeFigure();
        
//          Graph g = ea.matlabToGraph("graph.mat", "graph");
//          ea.graphToArff(g);
    }
}
