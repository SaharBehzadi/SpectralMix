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
public class WeightedMajorizationDisplClustered extends AbstractLayout<Integer, Integer> implements IterativeContext {

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
    double[][] distances;
    double[][] coord;
    double[][] weights;
    double[][] groundTruth; //ground truth coordinates
    int[][] clIds; //clusterIDs; numlevels x numObj
    int[] clId;
    int numCl;
    int[] numCls;
    boolean isoInit = false;
    boolean mdsInit = false;
    boolean groundTruthInit = false;
    boolean randomInit = true;
    int convergeCounter = 0;
    boolean varianceFree = false;
    Graph g;
    GraphCompression gk;
    // int d;
    int n;
    int m;
    Random r;
    //double factor; //how much to keep of old coordinate
    Dimension dd; //size of display
    private static double PRECISION = 1E-12;
    static int maxIteration = 1000;

    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    public void setXydata(Point2D[] xydata) {
        this.xydata = xydata;
        scaleToDisplaySize();
    }

    public WeightedMajorizationDisplClustered(Graph<Integer, Integer> g, int[] clId, int numCl) {
        super(g);
        this.g = g;
        this.clId = clId;
        this.numCl = numCl;
        r = new Random(20);
        //this.factor = factor;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;

    }

    public WeightedMajorizationDisplClustered(Graph<Integer, Integer> g, int[][] clId, int[] numCl) {
        super(g);
        this.g = g;
        this.clIds = clId;
        this.numCls = numCl;
        r = new Random(20);
        //this.factor = factor;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;

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

    public WeightedMajorizationDisplClustered(Graph<Integer, Integer> g, double factor) {
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
        currentIteration = 0;
        bestloop = -1;
        loop = 0;
        dd = getSize();
        currentIteration = 0;
        cheatCounter = 0;
        cheatVar = false;
        vertices = new int[n];
        for (int i = 0; i < n; i++) {
            vertices[i] = i;
        }

        coord = new double[n][d];
        //init coord with isomap
        if (isoInit) {

            distances = pathdist();
            coord = new Matrix(MDSJ.classicalScaling(distances)).transpose().getArrayCopy();
        }
        if (groundTruthInit) {
            coord = groundTruth;
            distances = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
                }
            }

        }

        if (randomInit) {
            Random r = new Random(1);
            for (int i = 0; i < coord.length; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = r.nextDouble();
                }
            }
            distances = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
                }
            }
        }

        if (mdsInit) {
            distances = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    if (g.isNeighbor(i, j)) {
                        distances[i][j] = distances[j][i] = 1.0;
                    } else {
                        distances[i][j] = distances[j][i] = 2.0;
                    }
                }
                coord = new Matrix(MDSJ.stressMinimization(distances)).transpose().getArrayCopy();
            }

        }

        weights = new double[n][n];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < i; j++) {
//                weights[i][j] = weights[j][i] = 1.0;
//            }
//        }

//        
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
        gk.setCoord(coord);
        loop = 0;
        bestCost = gk.mdlFunctionSimpleSigmoid();
        paramCost = gk.paramCostBic();

        System.out.println("pCostsBic: " + paramCost);
        savedBits = costWithoutEmbedding - (bestCost + paramCost);
        System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
        scaleToDisplaySize();
    }

    public void InitializeOld() {
        currentIteration = 0;
        bestloop = -1;
        loop = 0;
        dd = getSize();
        vertices = new int[n];
        for (int i = 0; i < n; i++) {
            vertices[i] = i;
        }
        coord = new double[n][d];
        for (int i = 0; i < n; i++) {
            coord[i][0] = r.nextDouble() * 15;
            coord[i][1] = r.nextDouble() * 15;
        }
        distances = pathdist();
        weights = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
//                if(pdist[i][j] == 1 )
//                    weights[i][j] = weights[j][i] = 1;
//                else
                //weights[i][j] = weights[j][i] = 1;
                //weights[i][j] = weights[j][i] = pdist[i][j] * pdist[i][j];
                weights[i][j] = weights[j][i] = 100.0 * Math.exp(-distances[i][j]);
                //weights[i][j] = weights[j][i] = 10.0 * Math.pow(pdist[i][j], -1);
            }
        }

        //majorization of pich: weights 1: 25,269
//           StressMinimization st = new StressMinimization(pdist, coord, weights);
//           StressMinimization.majorize(coord, pdist, weights, 100);
//           coord = st.getPositions();
        //isomap: mdsj - 24,592
//           double[][] bla = MDSJ.classicalScaling(pdist, 2);
//           coord = new Matrix(bla).transpose().getArrayCopy();
//           //stress minimization - das ist änlich wie isomap mit weights 1 - 24,652
        //initialized with the output of classical scaling.
        //weights 100: 24,155, weights 1000: 23,766, weights 100,000: 23,385, weights 1000000 -> 23,515
        double[][] bla = MDSJ.stressMinimization(distances, weights, 2);
        coord = new Matrix(bla).transpose().getArrayCopy();
        gk = new GraphCompression(g, coord);
        costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
//

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
        bestCost = gk.mdlFunction();
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

    private synchronized void updateOldSigmoid() {
        try {
            double[][] weights_new = new double[n][n];
            double sumWeights_new = 0.0;
            double sumWeights_old = 0.0;
            double[][] dist_new = new double[n][n];
            double minWeight = 1E-9;
            double minDist = 1E-6;

            //compute sigmoid
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                }
            }
            Sigmoid s = new Sigmoid(p);
//            s.mu = 1.5;
//          //  s.sigma = 0.3;
//            s.a = 1.0;
//            s.b = 0.0;
            //DEBUG
            if (cheatVar) {
                if (s.sigma < 1.0) {
                    s.sigma = 1.0;
                }
            }
            //DEBUG

            //compute new weights and dists by setting derivatives of MDS error function and cost function equal
            for (int i = 0; i < g.getVertexCount(); i++) {
                for (int j = 0; j < i; j++) {
                    double x = dist(coord[i], coord[j]);

                    //double aktW = ((2 * s.mu - 2 * x) * (s.a - s.b)) / (Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));
                    if (graph.isNeighbor(i, j)) {
//                        //(2 * (a - b) ^ 2) / (pi * sigma ^ 2 * exp((2 * (mu - x) ^ 2) / sigma ^ 2) * log(2) * (b + (a - b) * (erf((mu - x) / sigma) + 1) - 1) ^ 2) + ((2 * mu - 2 * x) * (a - b)) / (pi ^ (1 / 2) * sigma ^ 3 * exp((mu - x) ^ 2 / sigma ^ 2) * log(2) * (b + (a - b) * (erf((mu - x) / sigma) + 1) - 1)) // (2 * (a - b) ^ 2)/
//                        double z1 = 2.0 * Math.pow((s.a - s.b), 2);
//                        //pi*sigma^2*exp((2*(mu - x)^2)/sigma^2)*log(2)*(b + (a - b)*(erf((mu - x)/sigma) + 1) - 1)^2 
//                        double exponentN1 = 2 * Math.pow((s.mu - x), 2) / Math.pow(s.sigma, 2);
//                        double squareN1 = s.b + (s.a - s.b) * (Erf.erf((s.mu - x) / s.sigma) + 1) - 1;
//                        double n1 = Math.PI * Math.pow(s.sigma, 2) * Math.exp(exponentN1) * Math.log(2) * Math.pow(squareN1, 2);
//                        //((2 * mu - 2 * x) * (a - b)) 
//                        double z2 = ((2 * s.mu - 2 * x) * (s.a - s.b));
//                        //pi ^ (1 / 2) * sigma ^ 3 * exp((mu - x) ^ 2 / sigma ^ 2) * log(2) * (b + (a - b) * (erf((mu - x) / sigma) + 1) - 1)    
//                        double exponentN2 = Math.pow((s.mu - x), 2) / Math.pow(s.sigma, 2);
//                        double n2 = Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(exponentN2) * Math.log(2) * (s.b + (s.a - s.b) * (Erf.erf((s.mu - x) / s.sigma) + 1));
//                        double aktW = z1 / n1 + z2 / n2;

//                        weights_new[i][j] = weights_new[j][i] = aktW;
                        double wEDGE = (s.a - s.b) * Math.pow(Math.PI, -0.1e1 / 0.2e1) * (-x + s.mu) * Math.pow(s.sigma, -0.3e1) * Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1) * Math.sqrt(0.2e1) / ((s.a - s.b) * (0.1e1 / 0.2e1 + Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) / 0.2e1)
                                + s.b) / Math.log(0.2e1) / 0.4e1 + Math.pow(s.a - s.b, 0.2e1) / Math.PI * Math.pow(Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1), 0.2e1) * Math.pow(s.sigma, -0.2e1) * Math.pow((s.a - s.b) * (0.1e1
                                / 0.2e1 + Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) / 0.2e1) + s.b,
                                -0.2e1) / Math.log(0.2e1) / 0.4e1;
                        weights_new[i][j] = weights_new[j][i] = Math.max(wEDGE, minWeight);
//                        
//

                        if (weights_new[i][j] > minWeight) {
                            //-(2*(a - b))/(pi^(1/2)*sigma*exp((mu - x)^2/sigma^2)*log(2)*(b + (a - b)*(erf((mu - x)/sigma) + 1))*((4*(a - b)^2)/(pi*sigma^2*exp((2*(mu - x)^2)/sigma^2)*log(2)*(b + (a - b)*(erf((mu - x)/sigma) + 1))^2) + (2*(2*mu - 2*x)*(a - b))/(pi^(1/2)*sigma^3*exp((mu - x)^2/sigma^2)*log(2)*(b + (a - b)*(erf((mu - x)/sigma) + 1)))))
//                            double zz1 = -(2 * (s.a - s.b));
//                            double expp1 = Math.pow((s.mu - x), 2) / Math.pow(s.sigma, 2);
//                            double erfTerm = (s.mu - x) / s.sigma;
//                            //(pi*sigma^2*exp((2*(mu - x)^2)/sigma^2)*log(2)*(b + (a - b)*(erf((mu - x)/sigma) + 1))^2)
//                            double nnn1 = Math.PI * Math.pow(s.sigma, 2) * Math.exp(
//                                    //(pi^(1/2)*sigma*exp((mu - x)^2/sigma^2)*log(2)*(b + (a - b)*(erf((mu - x)/sigma) + 1))*((4*(a - b)^2)/(pi*sigma^2*exp((2*(mu - x)^2)/sigma^2)*log(2)*(b + (a - b)*(erf((mu - x)/sigma) + 1))^2) + (2*(2*mu - 2*x)*(a - b))/(pi^(1/2)*sigma^3*exp((mu - x)^2/sigma^2)*log(2)*(b + (a - b)*(erf((mu - x)/sigma) + 1))))
//                                    //double nn1 = Math.pow(Math.PI, 0.5) * s.sigma * Math.exp(expp1) * Math.log(2) * (s.b + (s.a - s.b) * Erf.erf(erfTerm) + 1) * (4 * Math.pow((s.a - s.b), 2));
//
//                                    double  edgeChange = (s.a - s.b) / (Math.pow(Math.PI, 0.5) * s.sigma * weights_new[i][j] * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));
//                            double dn = x + edgeChange;

                            double dEDGE = -(-0.2e1 * Math.log(0.2e1) * Math.sqrt(Math.PI)
                                    * Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) * s.a * s.sigma * wEDGE * x + 0.2e1 * Math.log(0.2e1) * Math.sqrt(Math.PI) * Erf.erf(Math.sqrt(0.2e1)
                                    * (-x + s.mu) / s.sigma / 0.2e1) * s.b * s.sigma * wEDGE * x - 0.2e1
                                    * Math.log(0.2e1) * Math.sqrt(Math.PI) * s.a * s.sigma * wEDGE * x - 0.2e1
                                    * Math.log(0.2e1) * Math.sqrt(Math.PI) * s.b * s.sigma * wEDGE * x + Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1)
                                    * Math.sqrt(0.2e1) * s.a - Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1) * Math.sqrt(0.2e1) * s.b) * Math.pow(Math.PI, -0.1e1 / 0.2e1) / s.sigma / (s.a * Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) - s.b * Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) + s.a + s.b) / Math.log(0.2e1) / wEDGE / 0.2e1;

                            dist_new[i][j] = dist_new[j][i] = dEDGE;
                            //dist_new[i][j] = dist_new[j][i] = Math.max(minDist, dEDGE);

                        } else {
                            dist_new[i][j] = dist_new[j][i] = x;
                        }

                    } else {
//                        //(2 * (a - b) ^ 2) / (pi * sigma ^ 2 * exp((2 * (mu - x) ^ 2) / sigma ^ 2) * log(2) * (b + (a - b) * (erf((mu - x) / sigma) + 1) - 1) ^ 2) + ((2 * mu - 2 * x) * (a - b)) / (pi ^ (1 / 2) * sigma ^ 3 * exp((mu - x) ^ 2 / sigma ^ 2) * log(2) * (b + (a - b) * (erf((mu - x) / sigma) + 1) - 1)) // (2 * (a - b) ^ 2)/
//                        double z1 = 2.0 * Math.pow((s.a - s.b), 2);
//                        //pi*sigma^2*exp((2*(mu - x)^2)/sigma^2)*log(2)*(b + (a - b)*(erf((mu - x)/sigma) + 1) - 1)^2 
//                        double exponentN1 = 2 * Math.pow((s.mu - x), 2) / Math.pow(s.sigma, 2);
//                        double squareN1 = s.b + (s.a - s.b) * (Erf.erf((s.mu - x) / s.sigma) + 1) - 1;
//                        double n1 = Math.PI * Math.pow(s.sigma, 2) * Math.exp(exponentN1) * Math.log(2) * Math.pow(squareN1, 2);
//                        //((2 * mu - 2 * x) * (a - b)) 
//                        double z2 = ((2 * s.mu - 2 * x) * (s.a - s.b));
//                        //pi ^ (1 / 2) * sigma ^ 3 * exp((mu - x) ^ 2 / sigma ^ 2) * log(2) * (b + (a - b) * (erf((mu - x) / sigma) + 1) - 1)    
//                        double exponentN2 = Math.pow((s.mu - x), 2) / Math.pow(s.sigma, 2);
//                        double n2 = Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(exponentN2) * Math.log(2) * (s.b + (s.a - s.b) * (Erf.erf((s.mu - x) / s.sigma) + 1) - 1);
//                        double aktW = z1 / n1 + z2 / n2;
//                        weights_new[i][j] = weights_new[j][i] = aktW;
//                        //weights_new[i][j] = weights_new[j][i] = aktW;

                        double wNO = -(s.a - s.b) * Math.pow(Math.PI, -0.1e1 / 0.2e1) * (-x + s.mu) * Math.pow(s.sigma, -0.3e1) * Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1) * Math.sqrt(0.2e1) / ((s.a - s.b) * (0.1e1 / 0.2e1 - Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) / 0.2e1)
                                + s.b) / Math.log(0.2e1) / 0.4e1 + Math.pow(s.a - s.b, 0.2e1) / Math.PI * Math.pow(Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1), 0.2e1) * Math.pow(s.sigma, -0.2e1) * Math.pow((s.a - s.b) * (0.1e1
                                / 0.2e1 - Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) / 0.2e1) + s.b,
                                -0.2e1) / Math.log(0.2e1) / 0.4e1;

                        weights_new[i][j] = weights_new[j][i] = Math.max(wNO, minWeight);
                        if (weights_new[i][j] > minWeight) {
//                            double edgeChange = (s.a - s.b) / (Math.pow(Math.PI, 0.5) * s.sigma * weights_new[i][j] * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));
//                            double dn = x - edgeChange;
//                            if (Double.isNaN(dn)) {
//                                System.out.println("m");
//                            }

//                        double dNO = -(-0.2e1 * Math.log(0.2e1) * Math.sqrt(Math.PI)
//                                * Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) * s.a * s.sigma * wNO * x
//                                + 0.2e1 * Math.log(0.2e1) * Math.sqrt(Math.PI) * Erf.erf(Math.sqrt(0.2e1)
//                                * (-x + s.mu) / s.sigma / 0.2e1) * s.b * s.sigma * wNO * x - 0.2e1
//                                * Math.log(0.2e1) * Math.sqrt(Math.PI) * s.a * s.sigma * wNO * x - 0.2e1
//                                * Math.log(0.2e1) * Math.sqrt(Math.PI) * s.b * s.sigma * wNO * x + Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1)
//                                * Math.sqrt(0.2e1) * s.a - Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1) * Math.sqrt(0.2e1) * s.b) * Math.pow(Math.PI, -0.1e1 / 0.2e1) / s.sigma / (s.a * Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) - s.b * Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) + s.a + s.b) / Math.log(0.2e1) / wNO / 0.2e1;
                            double dNO = -(-0.2e1 * Math.log(0.2e1) * Math.sqrt(Math.PI)
                                    * Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) * s.a * s.sigma * wNO * x
                                    + 0.2e1 * Math.log(0.2e1) * Math.sqrt(Math.PI) * Erf.erf(Math.sqrt(0.2e1)
                                    * (-x + s.mu) / s.sigma / 0.2e1) * s.b * s.sigma * wNO * x + 0.2e1
                                    * Math.log(0.2e1) * Math.sqrt(Math.PI) * s.a * s.sigma * wNO * x + 0.2e1
                                    * Math.log(0.2e1) * Math.sqrt(Math.PI) * s.b * s.sigma * wNO * x + Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1)
                                    * Math.sqrt(0.2e1) * s.a - Math.exp(-Math.pow(-x + s.mu, 0.2e1) * Math.pow(s.sigma, -0.2e1) / 0.2e1) * Math.sqrt(0.2e1) * s.b) * Math.pow(Math.PI, -0.1e1 / 0.2e1) / s.sigma / (s.a * Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) - s.b * Erf.erf(Math.sqrt(0.2e1) * (-x + s.mu) / s.sigma / 0.2e1) - s.a - s.b) / Math.log(0.2e1) / wNO / 0.2e1;

                            //dist_new[i][j] = dist_new[j][i] = Math.max(minDist, dNO);
                            dist_new[i][j] = dist_new[j][i] = dNO;

                        } else {
                            dist_new[i][j] = dist_new[j][i] = x;
                        }

                        //weight only links 
//                      weights_new[i][j] = weights_new[j][i] = 1.0;   
//                       dist_new[i][j] = dist_new[j][i] = x;
                    }
                    sumWeights_new += weights_new[i][j];
                    sumWeights_old += weights[i][j];

                }
            }
            //Test
            //  weights_new = scale(weights_new);
            //dist_new = scale(dist_new);
            DataUtils du = new DataUtils();

            String name = "results_" + loop;

//            String name1 = "dist_" + loop;
//            ea.writeDoubleToMatlab(dist_new, name1);
            //Test
            //do embedding
            StressMinimization sm = new StressMinimization(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
            String ss = sm.iterate(1);
            //System.out.println(ss);
            double[][] bla = sm.getPositions();

            //double[][] bla = MDSJ.stressMinimization(dist_new, weights_new, 2);
            double[][] coord_new = new Matrix(bla).transpose().getArrayCopy();
            gk.setCoord(coord_new);
            double aktCost = gk.mdlFunction();
            if (aktCost > lastCost) {
                if (cheatVar == false) {
                    cheatVar = true;
                } else {
                    cheatVar = false;
                }
            } else {
                cheatVar = false;
                cheatCounter = 0;
            }

            System.out.println(loop + ", " + aktCost);
            lastCost = aktCost;

            //copy coord, distance and weights
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = coord_new[i][j];
                }
            }

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    weights[i][j] = weights[j][i] = weights_new[i][j];
                    distances[i][j] = distances[j][i] = dist_new[i][j];
                }
            }

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    // bestDb[i][j] = coord[i][j];
                }
                xydata[i].setLocation(coord[i][0], coord[i][1]);
            }

            if (aktCost < bestCost) {
                bestCost = aktCost;
                savedBits = costWithoutEmbedding - (aktCost + paramCost);
                bestloop = loop;
                //System.out.println(loop + ", " + aktCost + " saved bits: " + savedBits);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < d; j++) {
                        bestDb[i][j] = coord[i][j];
                    }
                    //xydata[i].setLocation(coord[i][0], coord[i][1]);
                }
            }
            scaleToDisplaySize();
            //Test

//            String name2 = "coord_" + loop;
//            ea.writeDoubleToMatlab(coord, name2);
            //Test
            if (loop == 100) {
                String fName = name + ".mat";

                du.saveAsMatlab4(distances, weights, du.colScaleData(coord), du.colScaleData(groundTruth), "dist", "weights", "coord", "gt", fName);
            }

            currentIteration++;

            loop++;
        } catch (ConcurrentModificationException cme) {
        }

    }
    
      private synchronized void updateOldO() {
        try {
            if(loop%2 == 1 || loop%2 == 0){
            
            double[][] weights_new = new double[numCl][numCl];
//            double sumWeights_new = 0.0;
//            double sumWeights_old = 0.0;
            double[][] dist_new = new double[numCl][numCl];
            double minWeight = 1E-9;
            double minDist = 1E-6;

            //compute sigmoid
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                }
            }
            SimpleSigmoid s = new SimpleSigmoid(p);
//            if (s.sigma > 1) {
//                s.sigma = 1;

            double aktCost = s.costAllPairs(p);

            //System.out.println(loop + " " + s.mu + " " + s.sigma + " Cost=" + aktCost);
            // if(cheatVar)
            //  s.sigma = Math.max(s.sigma, 5.0);
            // System.out.println("variance before update: " + s.sigma);
            //compute new weights and dists by setting derivatives of MDS error function and cost function equal
            int correctionCounter = 0;

            //cluster centers
            double[][] centers = new double[numCl][2];
            int[] numObjCl = new int[numCl];
            for (int i = 0; i < g.getVertexCount(); i++) {
                centers[clId[i]][0] += coord[i][0];
                centers[clId[i]][1] += coord[i][1];
                numObjCl[clId[i]]++;
            }
            for (int i = 0; i < numCl; i++) {
                centers[i][0] /= numObjCl[i];
                centers[i][1] /= numObjCl[i];
            }

            for (int i = 0; i < g.getVertexCount(); i++) {
                for (int j = 0; j < i; j++) {
                    double x = dist(coord[i], coord[j]);

                    //double aktW = ((2 * s.mu - 2 * x) * (s.a - s.b)) / (Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));
                    double[] dw = s.parabola(x, graph.isNeighbor(i, j));
                    weights_new[clId[i]][clId[j]] = weights_new[clId[j]][clId[i]] += dw[1];
                    dist_new[clId[i]][clId[j]] = dist_new[clId[j]][clId[i]] += dw[1] * (dw[0] - x);
//                    if (weights_new[i][j] < minWeight) {
//                        weights_new[i][j] = weights_new[j][i] = minWeight;
//                        dist_new[i][j] = dist_new[j][i] = distances[i][j];
//                        correctionCounter++;
//                    }
//                    sumWeights_new += weights_new[i][j];
//                    sumWeights_old += weights[i][j];

                }
            }
            for (int i = 0; i < numCl; i++) {
                for (int j = 0; j < numCl; j++) {
                    dist_new[i][j] /= weights_new[i][j];
                    dist_new[i][j] += dist(centers[i], centers[j]);
                }
            }
            System.out.println(loop + " " + correctionCounter + " " + m + " " + s.sigma + " " + aktCost + "clustered");
            //Test
            //  weights_new = scale(weights_new);
            //dist_new = scale(dist_new);
            DataUtils du = new DataUtils();

            String name = "results_" + loop;

            String name1 = "dist_" + loop;
            String name2 = "weights_" + loop;
            IO ea = new IO();
            // ea.writeDoubleToMatlab(dist_new, name1);
            // ea.writeDoubleToMatlab(weights_new, name2);
            //Test

            //do embedding
//            StressMinimization sm = new StressMinimization(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
//            String ss = sm.iterate(1);
//            //System.out.println(ss);
//            double[][] bla = sm.getPositions();
            //own majorization
            WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dist_new, new Matrix(centers).transpose().getArrayCopy(), weights_new);
            wp.iterate();
            double[][] bla = wp.getPositions();

            //double[][] bla = MDSJ.stressMinimization(dist_new, weights_new, 2);
            double[][] centers_new = new Matrix(bla).transpose().getArrayCopy();
            for (int i = 0; i < numCl; i++) {
                centers_new[i][0] -= centers[i][0];
                centers_new[i][1] -= centers[i][1];
            }
            double[][] coord_new = new double[g.getVertexCount()][2];
            for (int i = 0; i < g.getVertexCount(); i++) {
                coord_new[i][0] = coord[i][0] + centers_new[clId[i]][0];
                coord_new[i][1] = coord[i][1] + centers_new[clId[i]][1];
            }
            gk.setCoord(coord_new);
            //double aktCost = gk.mdlFunctionSimpleSigmoid(false);

            String name3 = "coord_" + loop;
            // du.saveAsMatlab3(dist_new, weights_new, coord_new, "dist", "weight", "coord", name+".mat");
            // ea.writeDoubleToMatlab(coord, name3);

            //DEBUG
            //gk.setEdgeCosts();
//             Visualization v = new Visualization(g);
//             String title = "loop_" + Integer.toString(loop);
//             v.displayCoord(coord, gk.edgeCosts, title);
//            if (aktCost > lastCost) {
//           //     cheatVar = true;
//                
//                if (cheatVar == false) {
//                    cheatVar = true;
//                } else {
//                    cheatVar = false;
//                }
//                }
//             else {
//                cheatVar = false;
//                cheatCounter = 0;
//            }
//        
            //   System.out.println(loop + ", " + aktCost + " " + s.sigma);
            if (lastCost - aktCost < 1 && convergeCounter < 4) {
                convergeCounter++;
            }
            lastCost = aktCost;

            //copy coord, distance and weights
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = coord_new[i][j];
                }
            }

//            for (int i = 0; i < n; i++) {
//                for (int j = 0; j < i; j++) {
//                    weights[i][j] = weights[j][i] = weights_new[i][j];
//                    distances[i][j] = distances[j][i] = dist_new[i][j];
//                }
//            }

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    // bestDb[i][j] = coord[i][j];
                }
                xydata[i].setLocation(coord[i][0], coord[i][1]);
            }

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
            scaleToDisplaySize();
            //Test

//            String name2 = "coord_" + loop;
//            ea.writeDoubleToMatlab(coord, name2);
            //Test
//            if (loop == 1) {
//                String fName = name + ".mat";
//
//                du.saveAsMatlab4(distances, weights, du.colScaleData(coord), du.colScaleData(groundTruth), "dist", "weights", "coord", "gt", fName);
//            }
            currentIteration++;

            loop++;
            }
            else{
                 double[][] weights_new = new double[n][n];
            double sumWeights_new = 0.0;
            double sumWeights_old = 0.0;
            double[][] dist_new = new double[n][n];
            double minWeight = 1E-9;
            double minDist = 1E-6;


            //compute sigmoid
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                }
            }
            SimpleSigmoid s = new SimpleSigmoid(p);
//            if (s.sigma > 1) {
//                s.sigma = 1;


            double aktCost = s.costAllPairs(p);

            //System.out.println(loop + " " + s.mu + " " + s.sigma + " Cost=" + aktCost);
            // if(cheatVar)
            //  s.sigma = Math.max(s.sigma, 5.0);
            // System.out.println("variance before update: " + s.sigma);



            //compute new weights and dists by setting derivatives of MDS error function and cost function equal
            int correctionCounter = 0;

            for (int i = 0; i < g.getVertexCount(); i++) {
                for (int j = 0; j < i; j++) {
                    double x = dist(coord[i], coord[j]);


                    //double aktW = ((2 * s.mu - 2 * x) * (s.a - s.b)) / (Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));

                    double[] dw = s.parabola(x, graph.isNeighbor(i, j));
                    dist_new[i][j] = dist_new[j][i] = dw[0];
                    weights_new[i][j] = weights_new[j][i] = dw[1];
//                    if (weights_new[i][j] < minWeight) {
//                        weights_new[i][j] = weights_new[j][i] = minWeight;
//                        dist_new[i][j] = dist_new[j][i] = distances[i][j];
//                        correctionCounter++;
//                    }
                    sumWeights_new += weights_new[i][j];
                    sumWeights_old += weights[i][j];

                }
            }
            System.out.println(loop + " " + correctionCounter + " " + m + " " + s.sigma + " " + aktCost);
            //Test
            //  weights_new = scale(weights_new);
            //dist_new = scale(dist_new);
            DataUtils du = new DataUtils();

            String name = "results_" + loop;


            String name1 = "dist_" + loop;
            String name2 = "weights_" + loop;
            IO ea = new IO();
            // ea.writeDoubleToMatlab(dist_new, name1);
            // ea.writeDoubleToMatlab(weights_new, name2);
            //Test




            //do embedding
//            StressMinimization sm = new StressMinimization(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
//            String ss = sm.iterate(1);
//            //System.out.println(ss);
//            double[][] bla = sm.getPositions();

            //own majorization
            WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
            wp.iterate();
            double[][] bla = wp.getPositions();



            //double[][] bla = MDSJ.stressMinimization(dist_new, weights_new, 2);

            double[][] coord_new = new Matrix(bla).transpose().getArrayCopy();
            gk.setCoord(coord_new);
            //double aktCost = gk.mdlFunctionSimpleSigmoid(false);


            String name3 = "coord_" + loop;
            // du.saveAsMatlab3(dist_new, weights_new, coord_new, "dist", "weight", "coord", name+".mat");
            // ea.writeDoubleToMatlab(coord, name3);


            //DEBUG
            //gk.setEdgeCosts();
//             Visualization v = new Visualization(g);
//             String title = "loop_" + Integer.toString(loop);
//             v.displayCoord(coord, gk.edgeCosts, title);




//            if (aktCost > lastCost) {
//           //     cheatVar = true;
//                
//                if (cheatVar == false) {
//                    cheatVar = true;
//                } else {
//                    cheatVar = false;
//                }
//                }
//             else {
//                cheatVar = false;
//                cheatCounter = 0;
//            }
//        




            //   System.out.println(loop + ", " + aktCost + " " + s.sigma);
            if (lastCost - aktCost < 1 && convergeCounter < 4) {
                convergeCounter++;
            }
            lastCost = aktCost;

            //copy coord, distance and weights
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = coord_new[i][j];
                }
            }

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    weights[i][j] = weights[j][i] = weights_new[i][j];
                    distances[i][j] = distances[j][i] = dist_new[i][j];
                }
            }


            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    // bestDb[i][j] = coord[i][j];
                }
                xydata[i].setLocation(coord[i][0], coord[i][1]);
            }


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
            scaleToDisplaySize();
            //Test

//            String name2 = "coord_" + loop;
//            ea.writeDoubleToMatlab(coord, name2);
            //Test
//            if (loop == 1) {
//                String fName = name + ".mat";
//
//                du.saveAsMatlab4(distances, weights, du.colScaleData(coord), du.colScaleData(groundTruth), "dist", "weights", "coord", "gt", fName);
//            }


            currentIteration++;

            loop++;
            }
        } catch (ConcurrentModificationException cme) {
        }

    }

    private synchronized void update() {
        try {
            Random r = new Random();
            int index = r.nextInt(clIds.length);
            
            
            if(loop%2 == 1 || loop%2 == 0){
            
            double[][] weights_new = new double[numCls[index]][numCls[index]];
//            double sumWeights_new = 0.0;
//            double sumWeights_old = 0.0;
            double[][] dist_new = new double[numCls[index]][numCls[index]];
            double minWeight = 1E-9;
            double minDist = 1E-6;

            //compute sigmoid
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                }
            }
            SimpleSigmoid s = new SimpleSigmoid(p);
//            if (s.sigma > 1) {
//                s.sigma = 1;

            double aktCost = s.costAllPairs(p);

            //System.out.println(loop + " " + s.mu + " " + s.sigma + " Cost=" + aktCost);
            // if(cheatVar)
            //  s.sigma = Math.max(s.sigma, 5.0);
            // System.out.println("variance before update: " + s.sigma);
            //compute new weights and dists by setting derivatives of MDS error function and cost function equal
            int correctionCounter = 0;

            //cluster centers
            double[][] centers = new double[numCls[index]][2];
            int[] numObjCl = new int[numCls[index]];
            for (int i = 0; i < g.getVertexCount(); i++) {
                centers[clIds[index][i]][0] += coord[i][0];
                centers[clIds[index][i]][1] += coord[i][1];
                numObjCl[clIds[index][i]]++;
            }
            for (int i = 0; i < numCls[index]; i++) {
                centers[i][0] /= numObjCl[i];
                centers[i][1] /= numObjCl[i];
            }

            for (int i = 0; i < g.getVertexCount(); i++) {
                for (int j = 0; j < i; j++) {
                    double x = dist(coord[i], coord[j]);

                    //double aktW = ((2 * s.mu - 2 * x) * (s.a - s.b)) / (Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));
                    double[] dw = s.parabola(x, graph.isNeighbor(i, j));
                    weights_new[clIds[index][i]][clIds[index][j]] = weights_new[clIds[index][j]][clIds[index][i]] += dw[1];
                    dist_new[clIds[index][i]][clIds[index][j]] = dist_new[clIds[index][j]][clIds[index][i]] += dw[1] * (dw[0] - x);
//                    if (weights_new[i][j] < minWeight) {
//                        weights_new[i][j] = weights_new[j][i] = minWeight;
//                        dist_new[i][j] = dist_new[j][i] = distances[i][j];
//                        correctionCounter++;
//                    }
//                    sumWeights_new += weights_new[i][j];
//                    sumWeights_old += weights[i][j];

                }
            }
            for (int i = 0; i < numCls[index]; i++) {
                for (int j = 0; j < numCls[index]; j++) {
                    dist_new[i][j] /= weights_new[i][j];
                    dist_new[i][j] += dist(centers[i], centers[j]);
                }
            }
            System.out.println(loop + " " + correctionCounter + " " + m + " " + s.sigma + " " + aktCost + " " + index);
            //Test
            //  weights_new = scale(weights_new);
            //dist_new = scale(dist_new);
            DataUtils du = new DataUtils();

            String name = "results_" + loop;

            String name1 = "dist_" + loop;
            String name2 = "weights_" + loop;
            IO ea = new IO();
            // ea.writeDoubleToMatlab(dist_new, name1);
            // ea.writeDoubleToMatlab(weights_new, name2);
            //Test

            //do embedding
//            StressMinimization sm = new StressMinimization(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
//            String ss = sm.iterate(1);
//            //System.out.println(ss);
//            double[][] bla = sm.getPositions();
            //own majorization
            WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dist_new, new Matrix(centers).transpose().getArrayCopy(), weights_new);
            wp.iterate();
            double[][] bla = wp.getPositions();

            //double[][] bla = MDSJ.stressMinimization(dist_new, weights_new, 2);
            double[][] centers_new = new Matrix(bla).transpose().getArrayCopy();
            for (int i = 0; i < numCls[index]; i++) {
                centers_new[i][0] -= centers[i][0];
                centers_new[i][1] -= centers[i][1];
            }
            double[][] coord_new = new double[g.getVertexCount()][2];
            for (int i = 0; i < g.getVertexCount(); i++) {
                coord_new[i][0] = coord[i][0] + centers_new[clIds[index][i]][0];
                coord_new[i][1] = coord[i][1] + centers_new[clIds[index][i]][1];
            }
            gk.setCoord(coord_new);
            //double aktCost = gk.mdlFunctionSimpleSigmoid(false);

            String name3 = "coord_" + loop;
            // du.saveAsMatlab3(dist_new, weights_new, coord_new, "dist", "weight", "coord", name+".mat");
            // ea.writeDoubleToMatlab(coord, name3);

            //DEBUG
            //gk.setEdgeCosts();
//             Visualization v = new Visualization(g);
//             String title = "loop_" + Integer.toString(loop);
//             v.displayCoord(coord, gk.edgeCosts, title);
//            if (aktCost > lastCost) {
//           //     cheatVar = true;
//                
//                if (cheatVar == false) {
//                    cheatVar = true;
//                } else {
//                    cheatVar = false;
//                }
//                }
//             else {
//                cheatVar = false;
//                cheatCounter = 0;
//            }
//        
            //   System.out.println(loop + ", " + aktCost + " " + s.sigma);
            if (lastCost - aktCost < 1 && convergeCounter < 4) {
                convergeCounter++;
            }
            lastCost = aktCost;

            //copy coord, distance and weights
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = coord_new[i][j];
                }
            }

//            for (int i = 0; i < n; i++) {
//                for (int j = 0; j < i; j++) {
//                    weights[i][j] = weights[j][i] = weights_new[i][j];
//                    distances[i][j] = distances[j][i] = dist_new[i][j];
//                }
//            }

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    // bestDb[i][j] = coord[i][j];
                }
                xydata[i].setLocation(coord[i][0], coord[i][1]);
            }

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
            scaleToDisplaySize();
            //Test

//            String name2 = "coord_" + loop;
//            ea.writeDoubleToMatlab(coord, name2);
            //Test
//            if (loop == 1) {
//                String fName = name + ".mat";
//
//                du.saveAsMatlab4(distances, weights, du.colScaleData(coord), du.colScaleData(groundTruth), "dist", "weights", "coord", "gt", fName);
//            }
            currentIteration++;

            loop++;
            }
            else{
                 double[][] weights_new = new double[n][n];
            double sumWeights_new = 0.0;
            double sumWeights_old = 0.0;
            double[][] dist_new = new double[n][n];
            double minWeight = 1E-9;
            double minDist = 1E-6;


            //compute sigmoid
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                }
            }
            SimpleSigmoid s = new SimpleSigmoid(p);
//            if (s.sigma > 1) {
//                s.sigma = 1;


            double aktCost = s.costAllPairs(p);

            //System.out.println(loop + " " + s.mu + " " + s.sigma + " Cost=" + aktCost);
            // if(cheatVar)
            //  s.sigma = Math.max(s.sigma, 5.0);
            // System.out.println("variance before update: " + s.sigma);



            //compute new weights and dists by setting derivatives of MDS error function and cost function equal
            int correctionCounter = 0;

            for (int i = 0; i < g.getVertexCount(); i++) {
                for (int j = 0; j < i; j++) {
                    double x = dist(coord[i], coord[j]);


                    //double aktW = ((2 * s.mu - 2 * x) * (s.a - s.b)) / (Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));

                    double[] dw = s.parabola(x, graph.isNeighbor(i, j));
                    dist_new[i][j] = dist_new[j][i] = dw[0];
                    weights_new[i][j] = weights_new[j][i] = dw[1];
//                    if (weights_new[i][j] < minWeight) {
//                        weights_new[i][j] = weights_new[j][i] = minWeight;
//                        dist_new[i][j] = dist_new[j][i] = distances[i][j];
//                        correctionCounter++;
//                    }
                    sumWeights_new += weights_new[i][j];
                    sumWeights_old += weights[i][j];

                }
            }
            System.out.println(loop + " " + correctionCounter + " " + m + " " + s.sigma + " " + aktCost);
            //Test
            //  weights_new = scale(weights_new);
            //dist_new = scale(dist_new);
            DataUtils du = new DataUtils();

            String name = "results_" + loop;


            String name1 = "dist_" + loop;
            String name2 = "weights_" + loop;
            IO ea = new IO();
            // ea.writeDoubleToMatlab(dist_new, name1);
            // ea.writeDoubleToMatlab(weights_new, name2);
            //Test




            //do embedding
//            StressMinimization sm = new StressMinimization(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
//            String ss = sm.iterate(1);
//            //System.out.println(ss);
//            double[][] bla = sm.getPositions();

            //own majorization
            WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
            wp.iterate();
            double[][] bla = wp.getPositions();



            //double[][] bla = MDSJ.stressMinimization(dist_new, weights_new, 2);

            double[][] coord_new = new Matrix(bla).transpose().getArrayCopy();
            gk.setCoord(coord_new);
            //double aktCost = gk.mdlFunctionSimpleSigmoid(false);


            String name3 = "coord_" + loop;
            // du.saveAsMatlab3(dist_new, weights_new, coord_new, "dist", "weight", "coord", name+".mat");
            // ea.writeDoubleToMatlab(coord, name3);


            //DEBUG
            //gk.setEdgeCosts();
//             Visualization v = new Visualization(g);
//             String title = "loop_" + Integer.toString(loop);
//             v.displayCoord(coord, gk.edgeCosts, title);




//            if (aktCost > lastCost) {
//           //     cheatVar = true;
//                
//                if (cheatVar == false) {
//                    cheatVar = true;
//                } else {
//                    cheatVar = false;
//                }
//                }
//             else {
//                cheatVar = false;
//                cheatCounter = 0;
//            }
//        




            //   System.out.println(loop + ", " + aktCost + " " + s.sigma);
            if (lastCost - aktCost < 1 && convergeCounter < 4) {
                convergeCounter++;
            }
            lastCost = aktCost;

            //copy coord, distance and weights
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = coord_new[i][j];
                }
            }

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    weights[i][j] = weights[j][i] = weights_new[i][j];
                    distances[i][j] = distances[j][i] = dist_new[i][j];
                }
            }


            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    // bestDb[i][j] = coord[i][j];
                }
                xydata[i].setLocation(coord[i][0], coord[i][1]);
            }


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
            scaleToDisplaySize();
            //Test

//            String name2 = "coord_" + loop;
//            ea.writeDoubleToMatlab(coord, name2);
            //Test
//            if (loop == 1) {
//                String fName = name + ".mat";
//
//                du.saveAsMatlab4(distances, weights, du.colScaleData(coord), du.colScaleData(groundTruth), "dist", "weights", "coord", "gt", fName);
//            }


            currentIteration++;

            loop++;
            }
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

    private synchronized void updateOld() {
        try {
            double[][] weights_new = new double[n][n];
            double sumWeights_new = 0.0;
            double sumWeights_old = 0.0;
            Sigmoid s = new Sigmoid(0.5);

            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                }
            }
            //if (loop%10==0)
            s = new Sigmoid(p);
            double sumllh = 0;
            for (int i = 0; i < g.getVertexCount(); i++) {
                for (int j = 0; j < i; j++) {
                    double h = 0;
                    if (graph.isNeighbor(i, j)) {
                        h = -log2(s.f(dist(coord[i], coord[j])));
                    } else {
                        h = -log2(1.0 - s.f(dist(coord[i], coord[j])));
                    }
                    //weights[i][j] = weights[j][i] = weights[j][i] + h;
                    weights_new[i][j] = weights_new[j][i] = h;
                    sumWeights_new += h;
                    sumWeights_old += weights[i][j];
                    sumllh += h;
                }
            }
            double weightChange = sumWeights_new / sumWeights_old;
            double factor = Math.max(weightChange, maxWeightChange);
            double meanWeightsNew = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    weights[i][j] = weights[j][i] = weights[j][i] + factor * weights_new[i][j];
                    meanWeightsNew += (weights[i][j] / (double) m);

                }
            }
            double stdWeightsNew = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    stdWeightsNew += ((weights[i][j] - meanWeightsNew) * (weights[i][j] - meanWeightsNew)) / (double) m;
                }
            }
            stdWeightsNew = Math.sqrt(stdWeightsNew);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    double threshold = meanWeightsNew + stdWeightsNew;
                    if (weights[i][j] > threshold) {
                        weights[i][j] = weights[j][i] = 1.0;
                    }
                }
            }

            if (sumllh < bestCost) {
                bestCost = sumllh;
                savedBits = costWithoutEmbedding - (sumllh + paramCost);
                bestloop = loop;
                System.out.println(loop + ", " + sumllh + " saved bits: " + savedBits);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < d; j++) {
                        bestDb[i][j] = coord[i][j];
                    }
                    xydata[i].setLocation(coord[i][0], coord[i][1]);
                }
                scaleToDisplaySize();
                //Test
                IO ea = new IO();
                String name = "coord_" + loop;
                ea.writeDoubleToMatlab(coord, name);
                //Test

            }

            for (int innerloop = 0; innerloop < 10; innerloop++) {
                double[][] dbRestoredNew = new double[n][d];
                double[][] sij = new double[n][n];
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        sij[i][j] = 0;
                        double h = dist(coord[i], coord[j]);
                        if (h != 0) {
                            sij[i][j] = distances[i][j] / h;
                        }
                    }
                }
                for (int i = 0; i < n; i++) {
                    double sumWeight = 0;
                    for (int j = 0; j < n; j++) {
                        if (i != j) {
                            double h = weights[i][j]; //0;
                            //if (adjacency[i][j]) {
                            //    h = -log2(s.f(dist(dbRestored[i], dbRestored[j])));
                            //} else {
                            //    h = -log2(1.0 - s.f(dist(dbRestored[i], dbRestored[j])));
                            //}
                            sumWeight += h;
                            for (int jj = 0; jj < d; jj++) {
                                dbRestoredNew[i][jj] += h * (coord[j][jj] + sij[i][j] * (coord[i][jj] - coord[j][jj]));
                            }
                        }
                    }
                    for (int jj = 0; jj < d; jj++) {
                        dbRestoredNew[i][jj] /= sumWeight;
                    }
                }
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < d; j++) {
                        coord[i][j] = dbRestoredNew[i][j];
                    }
                }
            }

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
