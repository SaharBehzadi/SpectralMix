package nature;

//package edu.uci.ics.jung.algorithms.layout;
import Jama.Matrix;
import edu.uci.ics.jung.algorithms.layout.AbstractLayout;
import edu.uci.ics.jung.algorithms.layout.GraphElementAccessor;
import edu.uci.ics.jung.algorithms.layout.util.RandomLocationTransformer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
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
import java.util.List;
import java.util.Random;
import java.util.Vector;
import javax.swing.AbstractAction;
import javax.swing.JPopupMenu;
import mdsj.MDSJ;
import mdsj.StressMinimization;

/**
 * Implements the Kamada-Kawai algorithm for node layout.
 * Does not respect filter calls, and sometimes crashes when the view changes to it.
 *
 * @see "Tomihisa Kamada and Satoru Kawai: An algorithm for drawing general indirect graphs. Information Processing Letters 31(1):7-15, 1989"
 * @see "Tomihisa Kamada: On visualization of abstract objects and relations. Ph.D. dissertation, Dept. of Information Science, Univ. of Tokyo, Dec. 1988."
 *
 * @author Masanori Harada
 */
public class MultiModalMaj extends AbstractLayout<Integer, Integer> implements IterativeContext {

    private int currentIteration;
    private int maxTry = 100;
    private double maxWeightChange = 1.0;
    private int loop;
    private int bestloop;
    private double bestllh;
    private int d = 2;
    private String status = "OwnLayout";
    private boolean adjustForGravity = true;
    private int[] vertices;
    private Point2D[] xydata;
    double costWithoutEmbedding;
    double paramCost;
    double savedBits;
    Vector<Graph<Integer, Integer>> gg; //graphs
    double[][] dist; //Array of distances from the different modalities [numMod][m]
    double[][] numA; //Array of numerical coordinates in original space
    double[][] bestDb;
    double[][] pdist; //to be removed
    double[][] coord; //low dimensional coordinates
    double[][] weights; //Array of weights [numMod][m]
    double[][] jDist; //joint distance matrix
    double[][] jWeights; //weights to scale the joint distance matrix
    Graph g;
    GraphCompression gk;
    // int d;
    int n;
    int m;
    Random r;
    //double factor; //how much to keep of old coordinate
    Dimension dd; //size of display
    boolean verbose = true;

    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    public void setXydata(Point2D[] xydata) {
        this.xydata = xydata;
        scaleToDisplaySize();
    }

    public MultiModalMaj(Vector<Graph<Integer, Integer>> g, Vector<double[]> num, int visuIndex) {
        super(g.elementAt(visuIndex));
        this.g = g.elementAt(visuIndex);
        r = new Random(20);
        //this.factor = factor;
        n = g.elementAt(visuIndex).getVertexCount();
        m = n * (n - 1) / 2;
        dist = new double[g.size() + num.size()][m];
        weights = new double[g.size() + num.size()][m];
        numA = new double[num.size()][num.elementAt(0).length];
        jDist = new double[n][n];
        jWeights = new double[n][n];
        this.gg = g;
        //compute distances Graph: shortest path, numerical: Euclidean
        int counter = 0;
        for (int i = 0; i < g.size(); i++) {
            dist[counter] = pDist(g.elementAt(i));
            weights[counter] = pathWeights(100, counter, true);
            dist[counter] = scale(dist[counter]);
            counter++;
        }
        for (int i = 0; i < num.size(); i++) {
            dist[counter] = eDist(num.elementAt(i));
            dist[counter] = scale(dist[counter]);
            weights[counter] = ones();
            numA[i] = num.elementAt(i);
            center(i); //for numeric MDL
            counter++;
        }


        //check distance matrices
        if (verbose) {
            IO ea = new IO();
            ea.writeDoubleToMatlab(dist, "distanceCheck");
            ea.writeDoubleToMatlab(weights, "weightCheck");
        }

    }

    private void center(int index) {
        double mean = 0.0;
        for (int i = 0; i < numA[index].length; i++) {
            mean += (numA[index][i] / numA[index].length);
        }
        for (int i = 0; i < numA[index].length; i++) {
            numA[index][i] -= mean;
        }

    }

    private double[] ones() {
        double[] ones = new double[m];
        for (int i = 0; i < m; i++) {
            ones[i] = 0.5;
        }
        return ones;
    }

    private double[] pathWeights(double factor, int index, boolean scale) {
        double[] w = new double[m];
        double wMax = 0;
        for (int i = 0; i < m; i++) {
            w[i] = factor * Math.exp(-dist[index][i]);
            if (w[i] > wMax) {
                wMax = w[i];
            }
        }
        if (scale) {
            for (int i = 0; i < m; i++) {
                w[i] /= wMax;
            }
        }

        return w;
    }

    private double[] scale(double[] d) {
        double max = 0.0;
        for (int i = 0; i < d.length; i++) {
            if (d[i] > max) {
                max = d[i];
            }
        }
        for (int i = 0; i < d.length; i++) {
            d[i] /= max;
        }
        return d;
    }

    public double getBestllh() {
        return bestllh;
    }

    public double[][] getCoordinates() {
        double[][] c = new double[xydata.length][2];
        for (int i = 0; i < xydata.length; i++) {
            c[i][0] = xydata[i].getX();
            c[i][1] = xydata[i].getY();
        }
        return c;
    }

    public MultiModalMaj(Graph<Integer, Integer> g, double factor) {
        super(g);
        r = new Random(20);
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
        if (loop - bestloop < maxTry) {
            return false;
        } else {
            return true;
        }
    }

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }

    public void initialize() {
        currentIteration = 0;
        bestloop = -1;
        loop = 0;
        dd = getSize();
        currentIteration = 0;
        vertices = new int[n];
        for (int i = 0; i < n; i++) {
            vertices[i] = i;
        }

//        //embed each single distance matrix and determine cost
//        double[] costSingle = new double[dist.length];
//
//        for (int k = 0; k < dist.length; k++) {
//            double[][] distSquare = new double[n][n];
//            for (int i = 1; i < n; i++) {
//                for (int j = 0; j < i; j++) {
//                    distSquare[i][j] = distSquare[j][i] = dist[k][i * (i - 1) / 2 + j];
//                }
//            }
//                if (k < gg.size()) {
//                    double[][] weightsSquare = new double[n][n];
//                    double ww[] = pathWeights(100, k, false);
//                    for (int i = 1; i < n; i++) {
//                        for (int j = 0; j < i; j++) {
//                            weightsSquare[i][j] = weightsSquare[j][i] = ww[i * (i - 1) / 2 + j];
//
//                        }
//                    }
//                 double[][] coord = MDSJ.stressMinimization(distSquare, weightsSquare, 2);
//
//
//                }
// else {
//
// }
//            }


            //construct weighted distance matrix
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    double sumWeights = 0.0;
                    for (int k = 0; k < weights.length; k++) {
                        sumWeights += weights[k][i * (i - 1) / 2 + j];
                    }
                    double sumDist = 0.0;
                    for (int k = 0; k < dist.length; k++) {
                        sumDist += weights[k][i * (i - 1) / 2 + j] * dist[k][i * (i - 1) / 2 + j];
                    }
                    jDist[i][j] = jDist[j][i] = sumDist / sumWeights;
                    jWeights[i][j] = jWeights[j][i] = sumWeights;
                }
            }

            double[][] bla = MDSJ.stressMinimization(jDist, jWeights, 2);
            coord = new Matrix(bla).transpose().getArrayCopy();
            gk = new GraphCompression(g, coord);

            NumericCompression nc = new NumericCompression(numA, coord);
            double cc = nc.codingCost();
            costWithoutEmbedding = gk.codingCostNoEmbedding() + nc.costWithoutEmbedding();
            System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);



            xydata = new Point2D[n];
            bestDb = new double[n][d];
            for (int i = 0; i < n; i++) {
                xydata[i] = transform(vertices[i]);
                xydata[i].setLocation(coord[i][0], coord[i][1]);
                bestDb[i][0] = coord[i][0];
                bestDb[i][1] = coord[i][1];
            }
            gk.setCoord(coord);
            loop = 0;
            bestllh = gk.mdlFunction() + cc;
            paramCost = n * log2(gg.size() * n / 2 + numA.length);
            System.out.println("pCosts: " + paramCost);
            savedBits = costWithoutEmbedding - (bestllh + paramCost);
            System.out.println("init: coding Costs " + bestllh + " saved bits: " + savedBits);
            scaleToDisplaySize();




        }

    private

     double[][] pathdist() {
        DijkstraShortestPath<Integer, Integer> alg = new DijkstraShortestPath(g);
        double[][] dist = new double[n][n];
        double maxDist_mds = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i > j) {
                    List<Integer> l = alg.getPath(i, j);
                    if (l.size() > 0) {
                        dist[i][j] = l.size();
                        dist[j][i] = dist[i][j];
                        if (dist[i][j] > maxDist_mds) {
                            maxDist_mds = dist[i][j];
                        }
//                    } else {
//                        dist[i][j] = maxDist;
//                        dist[j][i] = maxDist;
//                    }
                    }
                }
            }
        }
        for (int i = 0; i < dist.length; i++) {
            for (int j = 0; j < dist.length; j++) {
                if (dist[i][j] == 0.0 && i != j) {
                    dist[i][j] = maxDist_mds + 1;
                }
            }
        }
        return dist;
    }

    private synchronized void update() {
        try {
            //determine new weights by costs
            int counter = 0;
            double sumllh = 0;
            double graphCost = 0.0;
            //in each iteration and each modality, change weights only by factor. Graphs----------
            for (int ii = 0; ii < gg.size(); ii++) {
                double[] weights_new = new double[m];
                double sumWeights_new = 0.0;
                double sumWeights_old = 0.0;
                Sigmoid s = new Sigmoid(0.5);

                PairOwn[] p = new PairOwn[m];
                for (int i = 1; i < n; i++) {
                    for (int j = 0; j < i; j++) {
                        p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), gg.elementAt(ii).isNeighbor(i, j));
                    }
                }
                //if (loop%10==0)
                s = new Sigmoid(p);

                for (int i = 1; i < n; i++) {
                    for (int j = 0; j < i; j++) {
                        double h = 0;
                        if (gg.elementAt(ii).isNeighbor(i, j)) {
                            h = -log2(s.f(dist(coord[i], coord[j])));
                        } else {
                            h = -log2(1.0 - s.f(dist(coord[i], coord[j])));
                        }
                        weights_new[i * (i - 1) / 2 + j] = h;
                        sumWeights_new += h;
                        sumWeights_old += weights[counter][i * (i - 1) / 2 + j];
                        sumllh += h;
                        graphCost += h;
                    }
                }
                double weightChange = sumWeights_new / sumWeights_old;
                double factor = Math.max(weightChange, maxWeightChange);

                for (int i = 1; i < n; i++) {
                    for (int j = 0; j < i; j++) {
                        weights[counter][i * (i - 1) / 2 + j] = weights[counter][i * (i - 1) / 2 + j] + factor * weights_new[i * (i - 1) / 2 + j];

                    }
                }
                counter++;
            }
//            //numerical attributes: one global cost
//            NumericCompression nc = new NumericCompression(numA, coord);
//            double cost = nc.codingCost();
//            sumllh += cost;
//            double sumWeights_old = weights[counter][0] * m;
//            double weightChange = cost / sumWeights_old;
//            double factor = Math.max(weightChange, maxWeightChange);
//            double nw = cost / m;
//            for (int ii = 0; ii < numA.length; ii++) {
//                for (int i = 1; i < n; i++) {
//                    for (int j = 0; j < i; j++) {
//                        weights[counter][i * (i - 1) / 2 + j] = weights[counter][i * (i - 1) / 2 + j] + factor * nw;
//
//                    }
//                }
//                counter++;
//            }
            //numerical attributes: individual weights
            double[] weights_new = new double[m];
            NumericCompression nc = new NumericCompression(numA, coord);
            double cost = nc.codingCost();
            sumllh += cost;
            double sumWeightsOld = 0.0;
            double sumWeightsNew = 0.0;
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    weights_new[i * (i - 1) / 2 + j] = nc.getCost(i, j);
                    sumWeightsOld += weights[counter][i * (i - 1) / 2 + j];
                    sumWeightsNew += weights_new[i * (i - 1) / 2 + j];

                }
            }
            double weightChange = sumWeightsNew / sumWeightsOld;
            double factor = Math.max(weightChange, maxWeightChange);
            for (int ii = 0; ii < numA.length; ii++) {
                for (int i = 1; i < n; i++) {
                    for (int j = 0; j < i; j++) {
                        weights[counter][i * (i - 1) / 2 + j] = weights[counter][i * (i - 1) / 2 + j] + factor * weights_new[i * (i - 1) / 2 + j];

                    }
                }
                counter++;
            }





            //check for improvement
            if (sumllh < bestllh) {
                bestllh = sumllh;
                savedBits = costWithoutEmbedding - (sumllh + paramCost);
                bestloop = loop;
                System.out.println(loop + ", " + sumllh + " graphCost " + graphCost + " numCost " + cost + " saved: " + savedBits);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < d; j++) {
                        bestDb[i][j] = coord[i][j];
                    }
                    xydata[i].setLocation(coord[i][0], coord[i][1]);
                }
                scaleToDisplaySize();
                //Test
                IO ea = new IO();
                String name = "jWeights_" + loop;
                ea.writeDoubleToMatlab(jWeights, name);
                //Test

            }
            //compute weighted similarity matrix and do MDS
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    double sumWeights = 0.0;
                    for (int k = 0; k < weights.length; k++) {
                        sumWeights += weights[k][i * (i - 1) / 2 + j];
                    }
                    double sumDist = 0.0;
                    for (int k = 0; k < dist.length; k++) {
                        sumDist += weights[k][i * (i - 1) / 2 + j] * dist[k][i * (i - 1) / 2 + j];
                    }
                    // jDist[i][j] = jDist[j][i] = sumDist / sumWeights;
                    jWeights[i][j] = jWeights[j][i] = sumWeights;
                }
            }

            double[][] bla = MDSJ.stressMinimization(jDist, jWeights, 2);
            coord = new Matrix(bla).transpose().getArrayCopy();


            loop++;
        } catch (ConcurrentModificationException cme) {
        }

    }

    public void step() {
        //randomShuffle();

        update();
        //TEST
//        String name = new Integer(currentIteration).toString();
//        IO ea = new IO();
//        ea.writeDoubleToMatlab(weights, "w_" + name);
//        currentIteration++;
        //TEST

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

    public void scaleToDisplaySizeLargest() {
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
        double max = Math.max(max_x, max_y);
        double min = Math.min(min_x, min_y);
        for (int i = 0; i < xydata.length; i++) {
            double x_new = ((xydata[i].getX() - min) / (max - min)) * 0.99 * width;
            double y_new = ((xydata[i].getY() - min) / (max_y - min)) * 0.99 * height;
            xydata[i].setLocation(x_new, y_new);
        }
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

    private double[] pDist(Graph<Integer, Integer> graph) {
        DijkstraShortestPath<Integer, Integer> alg = new DijkstraShortestPath(graph);
        double[] pd = new double[m];
        double maxPDist = 0.0;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                List<Integer> l = alg.getPath(i, j);
                double ddd = l.size();
                pd[i * (i - 1) / 2 + j] = ddd;
                if (ddd > maxPDist) {
                    maxPDist = ddd;
                }
            }
        }
        for (int i = 0; i < pd.length; i++) {
            if (pd[i] == 0) {
                pd[i] = maxPDist;
            }
        }
        return pd;
    }

    private double[] eDist(double[] d) {
        double[] ed = new double[m];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                ed[i * (i - 1) / 2 + j] = Math.sqrt(Math.pow((d[i] - d[j]), 2));
            }
        }
        return ed;
    }
}
