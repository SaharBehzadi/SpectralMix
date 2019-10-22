package nature;

//package edu.uci.ics.jung.algorithms.layout;
import edu.uci.ics.jung.algorithms.layout.AbstractLayout;
import edu.uci.ics.jung.algorithms.layout.GraphElementAccessor;
import edu.uci.ics.jung.algorithms.layout.util.RandomLocationTransformer;
import edu.uci.ics.jung.algorithms.shortestpath.Distance;
import edu.uci.ics.jung.algorithms.util.IterativeContext;


import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.util.ConcurrentModificationException;

//import edu.uci.ics.jung.algorithms.*; //import edu.uci.ics.jung.algorithms.GraphStatistics;
//import edu.uci.ics.jung.algorithms.IterativeContext;
//import edu.uci.ics.jung.algorithms.shortestpath.Distance;
//import edu.uci.ics.jung.algorithms.shortestpath.UnweightedShortestPath;
//import edu.uci.ics.jung.algorithms.util.RandomLocationTransformer;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.AbstractPopupGraphMousePlugin;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Random;
import javax.swing.AbstractAction;
import javax.swing.JPopupMenu;

/**
 * Implements the Kamada-Kawai algorithm for node layout.
 * Does not respect filter calls, and sometimes crashes when the view changes to it.
 *
 * @see "Tomihisa Kamada and Satoru Kawai: An algorithm for drawing general indirect graphs. Information Processing Letters 31(1):7-15, 1989"
 * @see "Tomihisa Kamada: On visualization of abstract objects and relations. Ph.D. dissertation, Dept. of Information Science, Univ. of Tokyo, Dec. 1988."
 *
 * @author Masanori Harada
 */
public class RandomMove extends AbstractLayout<Integer, Integer> implements IterativeContext {

    private int currentIteration;
    private int maxIterations = 100;
    private String status = "OwnLayout";
    private boolean adjustForGravity = true;
    private int[] vertices;
    private int[] vv; //vertices
    private Point2D[] init;
    private Point2D[] xydata;
    int d;
    int n;
    int m;
    Random r;
    double sig;
    double factor;
    Dimension dd;





    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    public void setXydata(Point2D[] xydata) {
        this.xydata = xydata;
        scaleToDisplaySize();
    }

    public RandomMove(Point2D[] init, Graph<Integer, Integer> g) {
        super(g);
        this.init = init;
         r = new Random(20);
        sig = 10;
        factor = 0.98;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;

    }



    public double[][] getCoordinates(){
        double[][] c = new double[xydata.length][2];
        for(int i = 0; i < xydata.length; i++){
            c[i][0] = xydata[i].getX();
            c[i][1] = xydata[i].getY();
        }
        return c;
    }


    public RandomMove(Graph<Integer, Integer> g) {
        super(g);
        r = new Random(20);
        sig = 0.05;
        factor = 0.98;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
    }

 

    /**

     *
     *
     *

    public String getStatus() {
    return status + this.getSize();
    }

    public void setMaxIterations(int maxIterations) {
    this.maxIterations = maxIterations;
    }

    /**
     * This one is an incremental visualization.
     */
    public boolean isIncremental() {
        return true;
    }

    /**
     * Returns true once the current iteration has passed the maximum count.
     */
    public boolean done() {
        if (currentIteration > maxIterations) {
            return true;
        }
        return false;
    }

    public void initialize() {
        currentIteration = 0;
        Graph<Integer, Integer> graph = getGraph();
        Dimension d = getSize();

        if (graph != null && d != null) {

            vertices = new int[n];
            for(int i = 0; i < n; i++)
                vertices[i] = i;

//            vertices = (Integer[]) graph.getVertices().toArray();
//            for(int i = 0; i < vertices.length; i++)
//                vertices[i].
            vv = new int[n];

            xydata = new Point2D[n];
            for (int i = 0; i < n; i++) {
                xydata[i] = transform(vertices[i]);
                xydata[i].setLocation(init[i].getX(), init[i].getY());
            }
            scaleToDisplaySize();

//            // assign IDs to all visible vertices
//            while (true) {
//                try {
//                    int index = 0;
//                    for (V v : graph.getVertices()) {
//                        Point2D xyd = transform(v);
//                        vertices[index] = v;
//                        xydata[index] = xyd;
//                        index++;
//                    }
//                    break;
//                } catch (ConcurrentModificationException cme) {
//                }
//            }
//        }
        }
    }

    private synchronized void update() {
        try {
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {

                    p[i * (i - 1) / 2 + j] = new PairOwn(xydata[i].distance(xydata[j]), graph.isNeighbor(i, j));
                }
            }
            Sigmoid s = new Sigmoid(p);
            for (int i = 0; i < n; i++) {
                // now determine individual cost of all edges and non-edges of node i
                Point2D newCoord = (Point2D) xydata[i].clone();
                double costBefore = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        if (graph.isNeighbor(vertices[i], vertices[j])) {
                            costBefore -= log2(s.f(newCoord.distance(xydata[j])));
                        } else {
                            costBefore -= log2(1.0 - s.f(newCoord.distance(xydata[j])));
                        }
                    }
                }
                double new_x = newCoord.getX() + r.nextGaussian()  * sig;
                double new_y = newCoord.getY() + r.nextGaussian()   * sig;
                newCoord.setLocation(new_x, new_y);
//            for (int j = 0; j < d; j++) {
//                newCoord[j] += r.nextGaussian() * sig;
//            }
                // now determine new cost
                double costAfter = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        if (graph.isNeighbor(i, j)) {
                            costAfter -= log2(s.f(newCoord.distance(xydata[j])));
                        } else {
                            costAfter -= log2(1.0 - s.f(newCoord.distance(xydata[j])));
                        }
                    }
                }
                if (costAfter < costBefore) {
                    xydata[i] = newCoord;
                }
            }
            if (sig < 0.00001) {
                s = s;
            }
            scaleToDisplaySize();
            currentIteration++;
            //System.out.println(currentIteration);
            sig *= factor;
        } catch (ConcurrentModificationException cme) {
        }

    }

    public void step() {
        //randomShuffle();
        update();
    }
//finally {
//		fireStateChanged();
//		}

    public static double log2(double x) {
        final double l2 = 1.0 / Math.log(2.0);
        return Math.log(x) * l2;
    }

    private synchronized void randomShuffle() {
        for (int i = 0; i < n; i++) {
            xydata[i].setLocation(r.nextDouble(), r.nextDouble());
        }
        scaleToDisplaySize();
    }

    /**
     * Shift all vertices so that the center of gravity is located at
     * the center of the screen.
     */
    public void scaleToDisplaySize() {
        dd = getSize();
        double height = dd.getHeight();
        double width = dd.getWidth();
        //scale coordinates between 0 and 1
        double max_x = -java.lang.Double.MAX_VALUE;
        double max_y = -java.lang.Double.MAX_VALUE;
        double min_x = -max_x;
        double min_y = -max_y;
        for (int i = 0; i < n; i++) {

            if (xydata[i].getX() > max_x) {
                max_x = xydata[i].getX();
            }
            if (xydata[i].getX() < min_x) {
                min_x = xydata[i].getX();
            }
            if (xydata[i].getY() > max_y) {
                max_y = xydata[i].getY();
            }
            if (xydata[i].getY() < min_y) {
                min_y = xydata[i].getY();
            }
        }

        for (int i = 0; i < xydata.length; i++) {
            double x_new = ((xydata[i].getX() - min_x) / (max_x - min_x)) * 0.99 * width;
            double y_new = ((xydata[i].getY() - min_y) / (max_y - min_y)) * 0.99 * height;
            xydata[i].setLocation(x_new, y_new);
        }
    }

    /* (non-Javadoc)
     * @see edu.uci.ics.jung.visualization.layout.AbstractLayout#setSize(java.awt.Dimension)
     */
    @Override
    public void setSize(Dimension size) {
        setInitializer(new RandomLocationTransformer<Integer>(size));
        super.setSize(size);
    }

    /**
     * Enable or disable gravity point adjusting.
     */
    public void setAdjustForGravity(boolean on) {
        adjustForGravity = on;
    }

    /**
     * Returns true if gravity point adjusting is enabled.
     */
    public boolean getAdjustForGravity() {
        return adjustForGravity;
    }

    public void reset() {
        currentIteration = 0;
    }
}
