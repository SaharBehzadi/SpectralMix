package nature;

//package edu.uci.ics.jung.algorithms.layout;
import Jama.Matrix;
import edu.uci.ics.jung.algorithms.layout.AbstractLayout;
import edu.uci.ics.jung.algorithms.layout.GraphElementAccessor;
import edu.uci.ics.jung.algorithms.layout.util.RandomLocationTransformer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.algorithms.shortestpath.Distance;
import edu.uci.ics.jung.algorithms.util.IterativeContext;
import org.apache.commons.math3.special.Erf;

import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.util.ConcurrentModificationException;

//import edu.uci.ics.jung.algorithms.*; //import edu.uci.ics.jung.algorithms.GraphStatistics;
//import edu.uci.ics.jung.algorithms.IterativeContext;
//import edu.uci.ics.jung.algorithms.shortestpath.Distance;
//import edu.uci.ics.jung.algorithms.shortestpath.UnweightedShortestPath;
//import edu.uci.ics.jung.algorithms.util.RandomLocationTransformer;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.AbstractPopupGraphMousePlugin;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.List;
import java.util.Random;
import javax.swing.AbstractAction;
import javax.swing.JPopupMenu;
import mdsj.MDSJ;
import mdsj.StressMinimization;

/**
 * Implements the Kamada-Kawai algorithm for node layout. Does not respect
 * filter calls, and sometimes crashes when the view changes to it.
 *
 * @see "Tomihisa Kamada and Satoru Kawai: An algorithm for drawing general
 * indirect graphs. Information Processing Letters 31(1):7-15, 1989"
 * @see "Tomihisa Kamada: On visualization of abstract objects and relations.
 * Ph.D. dissertation, Dept. of Information Science, Univ. of Tokyo, Dec. 1988."
 *
 * @author Masanori Harada
 */
public class IncrementalDispl extends AbstractLayout<Integer, Integer> implements IterativeContext {

    private int currentIteration;
    private int maxTry = 100;
    private double maxWeightChange = 5.0;
    private int loop;
    private int bestloop;
    private double bestCost; //overall best cost
    private double lastCost = Double.MAX_VALUE; //cost in previous interation
    private double lastSigma; // sigma in previous iteration.
    private boolean cheatVar; // enlarge variance in next iteration
    private int cheatCounter; //how often cheating can happen until converged
    private int d = 2;
    private String status = "OwnLayout";
    private boolean adjustForGravity = true;
    private int[] vertices;
    private Point2D[] xydata;
    double costWithoutEmbedding;
    double paramCost;
    double savedBits;
    double[][] bestDb;
//    double[][] distances;
    double[][] coord;
    // double[][] weights;
    double[][] groundTruth; //ground truth coordinates
    boolean randomInit = true;
    int convergeCounter = 0;
    boolean varianceFree = false;
    Graph g;
    GraphCompression gk;
    SimpleSigmoid s;
    Random r;
    // int d;
    int n;
    int m;
    Dimension dd; //size of display
    //new stuff
    double[][] nomEdge;
    double[][] nomEdgeR;
    double[] denomEdge;
    int[] vNonEdge1;
    int[] vNonEdge2;
    double[][] nomNonEdge1;
    double[][] nomNonEdge1R;
    double[] denomNonEdge1;
    double[][] nomNonEdge2;
    double[][] nomNonEdge2R;
    double[] denomNonEdge2;
    double[][] nomVertex;
    double[] denomVertex;
    boolean[] valid;
//new
    private static double PRECISION = 1E-12;
    static int maxIteration = 5000;

    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    public void setXydata(Point2D[] xydata) {
        this.xydata = xydata;
        scaleToDisplaySize();
    }

    public IncrementalDispl(Graph<Integer, Integer> g) {
        super(g);
        this.g = g;
        r = new Random(20);
        //this.factor = factor;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        valid = new boolean[g.getEdgeCount()];
        nomEdge = new double[g.getEdgeCount()][d];
        nomEdgeR = new double[g.getEdgeCount()][d];
        denomEdge = new double[g.getEdgeCount()];
        vNonEdge1 = new int[g.getEdgeCount()];
        vNonEdge2 = new int[g.getEdgeCount()];
        nomNonEdge1 = new double[g.getEdgeCount()][d];
        nomNonEdge1R = new double[g.getEdgeCount()][d];
        denomNonEdge1 = new double[g.getEdgeCount()];
        nomNonEdge2 = new double[g.getEdgeCount()][d];
        nomNonEdge2R = new double[g.getEdgeCount()][d];
        denomNonEdge2 = new double[g.getEdgeCount()];
        nomVertex = new double[n][d];
        denomVertex = new double[n];

    }

    public double getBestllh() {
        return bestCost;
    }

    public double[][] getCoordinates() {

        return bestDb;
    }

    public void setGroundTruth(double[][] groundTruth) {
        this.groundTruth = groundTruth;
    }

    public IncrementalDispl(Graph<Integer, Integer> g, double factor) {
        super(g);
        r = new Random(20);
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
    }

    /**
     *
     *
     *
     *
     *
     * public String getStatus() { return status + this.getSize(); }
     *
     * public void setMaxIterations(int maxIterations) { this.maxIterations =
     * maxIterations; }
     *
     * /**
     * This one is an incremental visualization.
     */
    public boolean isIncremental() {
        return true;
    }

    /**
     * Returns true once the current iteration has passed the maximum count.
     */
//    public boolean done() {
//        if ((loop - bestloop) < maxTry) {
//            return false;
//        } else {
//            return true;
//        }
//    }
//    
    public boolean done() {
        if (currentIteration > maxIteration) {
            return true;
        }
        return false;
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
        // r = new Random(1);
        currentIteration = 0;
        bestloop = -1;
        loop = 0;
        dd = getSize();
        cheatCounter = 0;
        cheatVar = false;
        vertices = new int[n];
        for (int i = 0; i < n; i++) {
            vertices[i] = i;
        }

        coord = new double[n][d];

        if (randomInit) {
            Random r = new Random(2);
            for (int i = 0; i < coord.length; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = r.nextDouble();
                }
            }
//            DataUtils du = new DataUtils();
//            du.saveAsMatlab(coord, "randomInit", "randomInit.mat");

        }

        gk = new GraphCompression(g, coord);
        costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);

        xydata = new Point2D[n];
        bestDb = new double[n][d];
        for (int i = 0; i < n; i++) {
            xydata[i] = transform(vertices[i]);
            xydata[i].setLocation(coord[i][0], coord[i][1]);
            bestDb[i][0] = coord[i][0];
            bestDb[i][1] = coord[i][1];
        }
        computeSigmoid();
        gk.setCoord(coord);
        loop = 0;
        bestCost = gk.mdlFunctionSimpleSigmoid();
        paramCost = gk.paramCostBic();

        System.out.println("pCostsBic: " + paramCost);
        savedBits = costWithoutEmbedding - (bestCost + paramCost);
        System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
        scaleToDisplaySize();
    }

    private double[][] pathdist() {
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

    //scales symmetric distance matrix and weight matrix between 0...1. Assumes positive values only
    private double[][] scale(double[][] d) {
        double[][] ds = new double[n][n];
        double max = -Double.MAX_VALUE;
        double min = Double.MAX_VALUE;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (d[i][j] > max) {
                    max = d[i][j];
                }
                if (d[i][j] < min) {
                    min = d[i][j];
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double s = (d[i][j] - min) / (max - min);
                s *= 100;
                ds[i][j] = ds[j][i] = Math.max(s, 0.000000001);
            }
        }

        return ds;
    }

    private void computeSigmoid() {
        //compute sigmoid
        PairOwn[] p = new PairOwn[m];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
            }
        }
        s = new SimpleSigmoid(p);
    }

    private double getAvgDist() {
        double avgDist = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i < j) {
                    avgDist += dist(coord[i], coord[j]) / n;
                }
            }
        }
        return avgDist;
    }

    private double[] zij(double dij, double wij, int i, int j) {
        double[] res = new double[d];
        double sij = 0.0;
        for (int k = 0; k < d; k++) {
            sij += Math.pow((coord[i][k] - coord[j][k]), 2);
        }
        if (sij != 0) {
            sij = dij / Math.sqrt(sij);
        }
        for (int k = 0; k < d; k++) {
            res[k] = wij * (coord[j][k] + sij * (coord[i][k] - coord[j][k]));
        }
        return res;
    }

    private synchronized void update() {
        try {
//           //compute sigmoid
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                }
            }
             if (currentIteration % 100 == 0) {
            computeSigmoid();
            System.out.println(currentIteration + " " + s.sigma);
            }
            s.sigma = Math.min(s.sigma, 1.0);
            double aktCost = s.costAllPairs(p);
            //System.out.println(currentIteration + " " + s.sigma + " " + aktCost);
            r = new Random(currentIteration);

            for (int i = 0; i < g.getEdgeCount(); i++) {
                Pair<Integer> akt = g.getEndpoints(i);
                //System.out.println(akt.getFirst() + " " + akt.getSecond());
                int notEdge1 = r.nextInt(n);
                while (akt.getFirst() == notEdge1 || g.isNeighbor(akt.getFirst(), notEdge1)) {
                    notEdge1 = r.nextInt(n);
                }
                int notEdge2 = r.nextInt(n);
                while (akt.getSecond() == notEdge2 || g.isNeighbor(akt.getSecond(), notEdge2)) {
                    notEdge2 = r.nextInt(n);
                }
                if (currentIteration > 0) {
                    for (int j = 0; j < d; j++) {
                        coord[akt.getFirst()][j] = nomVertex[akt.getFirst()][j] / denomVertex[akt.getFirst()];
                        coord[akt.getSecond()][j] = nomVertex[akt.getSecond()][j] / denomVertex[akt.getSecond()];
                        coord[notEdge1][j] = nomVertex[notEdge1][j] / denomVertex[notEdge1];
                        coord[notEdge2][j] = nomVertex[notEdge2][j] / denomVertex[notEdge2];
                    }
                }
                double x = dist(coord[akt.getFirst()], coord[akt.getSecond()]);
                double[] dw = s.parabola(x, true);
                double[] h = zij(dw[0], dw[1], akt.getFirst(), akt.getSecond());
                for (int j = 0; j < d; j++) {
                    nomVertex[akt.getFirst()][j] += h[j] - (valid[i] ? nomEdge[i][j] : 0);
                }
                nomEdge[i] = h;
                denomVertex[akt.getFirst()] += dw[1] - (valid[i] ? denomEdge[i] : 0);
                h = zij(dw[0], dw[1], akt.getSecond(), akt.getFirst());
                for (int j = 0; j < d; j++) {
                    nomVertex[akt.getSecond()][j] += h[j] - (valid[i] ? nomEdgeR[i][j] : 0);
                }
                nomEdgeR[i] = h;
                denomVertex[akt.getSecond()] += dw[1] - (valid[i] ? denomEdge[i] : 0);
                denomEdge[i] = dw[1];

                x = dist(coord[akt.getFirst()], coord[notEdge1]);
                dw = s.parabola(x, false);
                h = zij(dw[0], dw[1], akt.getFirst(), notEdge1);
                for (int j = 0; j < d; j++) {
                    nomVertex[akt.getFirst()][j] += h[j] - (valid[i] ? nomNonEdge1[i][j] : 0);
                }
                nomNonEdge1[i] = h;
                denomVertex[akt.getFirst()] += dw[1] - (valid[i] ? denomNonEdge1[i] : 0);
                h = zij(dw[0], dw[1], notEdge1, akt.getFirst());
                for (int j = 0; j < d; j++) {
                    nomVertex[notEdge1][j] += h[j];// - (valid[i] ? nomNonEdge1R[i][j] : 0);
                    if(valid[i])
                        nomVertex[vNonEdge1[i]][j] -= nomNonEdge1R[i][j];
                }
                nomNonEdge1R[i] = h;
                denomVertex[notEdge1] += dw[1];// - (valid[i] ? denomNonEdge1[i] : 0);
                if(valid[i])
                    denomVertex[vNonEdge1[i]] -= denomNonEdge1[i];
                denomNonEdge1[i] = dw[1];
                vNonEdge1[i] = notEdge1;

                x = dist(coord[akt.getSecond()], coord[notEdge2]);
                //double aktW = ((2 * s.mu - 2 * x) * (s.a - s.b)) / (Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));
                dw = s.parabola(x, false);
                h = zij(dw[0], dw[1], akt.getSecond(), notEdge2);
                for (int j = 0; j < d; j++) {
                    nomVertex[akt.getSecond()][j] += h[j] - (valid[i] ? nomNonEdge2[i][j] : 0);
                }
                nomNonEdge2[i] = h;
                denomVertex[akt.getSecond()] += dw[1] - (valid[i] ? denomNonEdge2[i] : 0);
                h = zij(dw[0], dw[1], notEdge2, akt.getSecond());
                for (int j = 0; j < d; j++) {
                    nomVertex[notEdge2][j] += h[j] ;
                    if(valid[i])
                        nomVertex[vNonEdge2[i]][j] -= nomNonEdge2R[i][j];
                }
                nomNonEdge2R[i] = h;
                denomVertex[notEdge2] += dw[1];
                if(valid[i])
                    denomVertex[vNonEdge2[i]] -= denomNonEdge2[i];
                denomNonEdge2[i] = dw[1];
                vNonEdge2[i] = notEdge2;

                valid[i] = true;

            }

            DataUtils du = new DataUtils();
           
            if (lastCost - aktCost < 1 && convergeCounter < 4) {
                convergeCounter++;
            }
            lastCost = aktCost;
            //du.saveAsMatlab(coord, name3, "bla.mat");
            for (int i = 0; i < n; i++) {
                xydata[i].setLocation(coord[i][0], coord[i][1]);
            }
            scaleToDisplaySize();
            if (aktCost < bestCost) {
                bestCost = aktCost;
                savedBits = costWithoutEmbedding - (aktCost + paramCost);
                bestloop = loop;
                //System.out.println(loop + ", " + aktCost + " saved bits: " + savedBits);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < d; j++) {
                        bestDb[i][j] = coord[i][j];
                    }
                    //xydata[i].setLocation(coord[i][0], coord[i]e[1]);
                }
            }

            //Test
            String name22 = "coord_" + loop;
            //  ea.writeDoubleToMatlab(coord, name22);
 
            currentIteration++;

            loop++;
        } catch (ConcurrentModificationException cme) {
        }

    }

    public double getDistThreshold(double sigma) {
        double w = 1e-6;
        double arg = -Math.PI * Math.pow(Math.log(2), 2) * Math.pow(sigma, 4) * Math.pow(w, 2);
        double lower = Math.exp(0.5 * LambertW(arg));
        double upper = Math.sqrt(2) * w * Math.pow(sigma, 3) * Math.sqrt(Math.PI) * Math.log(2) - 1.5 * lower;
        return -upper / lower;

    }

    public static double LambertW(double z) {
        double S = 0.0;
        for (int n = 1; n <= 100; n++) {
            double Se = S * Math.pow(Math.E, S);
            double S1e = (S + 1)
                    * Math.pow(Math.E, S);
            if (PRECISION > Math.abs((z - Se) / S1e)) {
                return S;
            }
            S
                    -= (Se - z) / (S1e - (S + 2) * (Se - z) / (2 * S + 2));
        }
        return S;
    }

    public void step() {
        //randomShuffle();
//        while (!done()) {
//           
//        }
        //TEST
        update();

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
     * Shift all vertices so that the center of gravity is located at the center
     * of the screen.
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
            double x_new = ((xydata[i].getX() - min_x) / (max_x - min_x)) * 0.99 * width + 20;
            double y_new = ((xydata[i].getY() - min_y) / (max_y - min_y)) * 0.99 * height + 20;
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
