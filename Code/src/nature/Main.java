/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.decorators.EllipseVertexShapeTransformer;
import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.util.Vector;
import javax.swing.JFrame;
import org.apache.commons.collections15.Transformer;


/**
 *
 * @author claudia
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        IO ea = new IO();
        
        Graph g = ea.matlabToGraph("graph.mat", "graph");
        double[][] xx =  ea.readMatlabMatrix("coord.mat", "coord");
      double[][] ids = ea.readMatlabMatrix("id.mat", "id");
        int[] id = new int[ids.length];
        for (int i=0; i<ids.length ; i++)
            id[i] = (int) ids[i][0] ;
        
//        double[] x = ea.readMatlabMatrix("matlab.mat", "test")[0];
//        double[] y = ea.readMatlabMatrix("matlab.mat", "test")[1];
//        double[] z =  ea.readMatlabMatrix("matlab.mat", "test")[2];
//                
//
////
//        Vector<Graph<Integer, Integer>> gg = new  Vector<Graph<Integer, Integer>>();
//        gg.add(g);
//        Vector<double[]> num = new Vector<double[]>();
//        num.add(x);
//        num.add(y);
        //num.add(z);

         //Graph g = ea.matlabToGraph("adj.mat", "adj");

        Visualization v = new Visualization(g);
     double[][] coord = v.getCoordinatesIsomapOnly();
        Matrix m = new Matrix(xx).transpose();
        double[][] gs = m.getArrayCopy();
        v.displayCoordSmall(coord,"Gold Standard", id);
        GraphCompression gk = new GraphCompression(g, xx);
        System.out.println(gk.codingCostNoEmbedding());
        System.out.println(gk.mdlFunction());
         //v.getCoordinatesItMaj();
        //v.getCoordinatesMultiModalMaj(gg, num, 0);




       

    

//        Summarization s = new Summarization(g);
//        s.graphToDistanceMatrix();
//        s.greedyMerge();

        //multi-threaded
//        int[] ids = new int[g.getVertexCount()];
//        String[] s = new String[g.getVertexCount()];
//        Visualization v = new Visualization(g);
//        v.getCoordinatesOwn(ids, s);

        //v.getCoordinatesOwn();

//
//        Embedding e = new Embedding(g, 2);
//        //e.setCoord(new Matrix(c).transpose().getArrayCopy());
//        e.randomInit(17);
//        e.isoInit();
////        e.ownInit();
//
//        // e.distanceScaling(5000);
//        e.improveCoordinates();


    }
}
