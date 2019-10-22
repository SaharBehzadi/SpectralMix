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
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.Vector;
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
public class WeightedMajorizationDisplGrid extends AbstractLayout<Integer, Integer> implements IterativeContext {

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
    double mu; //store parameters from the previous iteration in order to determine threshold distance
    double sigma;
    double tau; //threshold for weights to be considered
    double[][] bestDb;
    double[][] distances;
    double[][] positions;
    double[][] coord;
    double[][] weights;
    double[][] groundTruth; //ground truth coordinates
    boolean isoInit = false;
    boolean mdsInit = false;
    boolean groundTruthInit = false;
    boolean randomInit = true;
    int convergeCounter = 0;
    boolean varianceFree = true;
    Graph g;
    GraphCompression gk;
    // int d;
    int n;
    int m;
    Random r;
    //double factor; //how much to keep of old coordinate
    Dimension dd; //size of display
    int numBins;
    private static double PRECISION = 1E-12;
    static int maxIteration = 5000;

    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    public void setXydata(Point2D[] xydata) {
        this.xydata = xydata;
        scaleToDisplaySize();
    }

    public WeightedMajorizationDisplGrid(Graph<Integer, Integer> g) {
        super(g);
        this.g = g;
        r = new Random(21);
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

    public WeightedMajorizationDisplGrid(Graph<Integer, Integer> g, double factor) {
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
        positions = new double[n][d];
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
            Random r = new Random(10);
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
        mu = gk.s.mu;
        sigma = gk.s.sigma;
        tau = 0.0001;
        System.out.println("pCostsBic: " + paramCost);
        savedBits = costWithoutEmbedding - (bestCost + paramCost);
        System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
        scaleToDisplaySize();
    }

    public void initialize1603() {
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
        positions = new double[n][d];
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
            Random r = new Random(10);
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

        DataUtils du = new DataUtils();
        coord = du.colScaleData(coord);
        // du.saveAsMatlab(coord, "coord", "scaled.mat"); -- checked ok


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
        mu = gk.s.mu;
        sigma = gk.s.sigma;
        tau = 1E-9;

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

    private synchronized void update() {
        try {
            double[][] weights_new = new double[n][n];
            double sumWeights_new = 0.0;
            double sumWeights_old = 0.0;
            double[][] dist_new = new double[n][n];
            double minWeight = 1E-6;
            double minDist = 1E-6;
            DataUtils du = new DataUtils();


            double cutOff = gridThreshold();
            //cutOff = 100;
            double[][] minMax = du.minMax(coord);
            int numBinsX = (int) Math.ceil((minMax[0][1] - minMax[0][0]) / cutOff);
            int numBinsY = (int) Math.ceil((minMax[1][1] - minMax[1][0]) / cutOff);
           System.out.println(cutOff + " " + numBinsX + " " + numBinsY);


            int[][] numPoints = new int[numBinsX][numBinsY];
            for (int i = 0; i < n; i++) {
                //compute grid cell of point
                //first bin: [0, 0.01[, last bin [0.9, 1.0]
                int iIndex = Math.min((int) Math.floor((coord[i][0] - minMax[0][0]) / cutOff), numBinsX - 1);
                int jIndex = Math.min((int) Math.floor((coord[i][1] - minMax[1][0]) / cutOff), numBinsY - 1);
//              
                numPoints[iIndex][jIndex]++;
            }

            int[][][] pId = new int[numBinsX][numBinsY][];
            for (int i = 0; i < numBinsX; i++) {
                for (int j = 0; j < numBinsY; j++) {
                    pId[i][j] = new int[numPoints[i][j]];
                }
            }


            for (int i = 0; i < numBinsX; i++) {
                for (int j = 0; j < numBinsY; j++) {
                    numPoints[i][j] = 0;
                }
            }

            for (int i = 0; i < n; i++) {
                //compute grid cell of point
                //first bin: [0, 0.01[, last bin [0.9, 1.0]
                int iIndex = Math.min((int) Math.floor((coord[i][0] - minMax[0][0]) / cutOff), numBinsX - 1);
                int jIndex = Math.min((int) Math.floor((coord[i][1] - minMax[1][0]) / cutOff), numBinsY - 1);
//              
                pId[iIndex][jIndex][numPoints[iIndex][jIndex]] = i;
                numPoints[iIndex][jIndex]++;
            }
          
//
//            //store the point id in each grid cell
//            Vector<Integer>[][] points = (Vector<Integer>[][]) new Vector[numBinsX][numBinsY];
//            for (int i = 0; i < points.length; i++) {
//                for (int j = 0; j < points[i].length; j++) {
//                    points[i][j] = new Vector<Integer>();
//                }
//            }
//
//
//            for (int i = 0; i < n; i++) {
//                //compute grid cell of point
//                //first bin: [0, 0.01[, last bin [0.9, 1.0]
//                int iIndex = Math.min((int) Math.floor((coord[i][0] - minMax[0][0]) / cutOff), numBinsX - 1);
//                int jIndex = Math.min((int) Math.floor((coord[i][1] - minMax[1][0]) / cutOff), numBinsY - 1);
////              
//                points[iIndex][jIndex].add(i);
//            }
//
//



            //compute sigmoid
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                    //DEBUG
                    if (Double.isNaN(p[i * (i - 1) / 2 + j].dist)) {
                        //    System.out.println("dist NaN");
                    }
                    //DEBUG
                }
            }

            SimpleSigmoid s = new SimpleSigmoid(p);
            mu = s.mu;
            sigma = s.sigma;



            double aktCost = s.costAllPairs(p);

            //System.out.println(loop + " " + s.mu + " " + s.sigma + " Cost=" + aktCost);
            // if(cheatVar)
            //  s.sigma = Math.max(s.sigma, 5.0);
            // System.out.println("variance before update: " + s.sigma);



            //compute new weights and dists by setting derivatives of MDS error function and cost function equal
            int correctionCounter = 0;

            //double[] sumWeights = new double[n];
            //jeden punkt einzeln

            double[] sumWeights = new double[n];
            
//            //DEBUG
//            int[][] dwUpdated = new int[n][n];
//            //DEBUG
            
            //DEBUG
            
            //DEBUG
            
            positions = new double[n][d];

            //objects in the same grid cell: pId[i][j][b]
            for (int i = 0; i < numBinsX; i++) {
                for (int j = 0; j < numBinsY; j++) {
                    for (int a = 0; a < pId[i][j].length; a++) {
                        for (int b = 0; b < a; b++) {
                            double x = dist(coord[pId[i][j][a]], coord[pId[i][j][b]]);
                            if (!g.isNeighbor(pId[i][j][a], pId[i][j][b])) {
//                                dwUpdated[pId[i][j][a]][pId[i][j][b]]++;
//                                dwUpdated[pId[i][j][b]][pId[i][j][a]]++;
                                double[] dw = s.parabola(x, false);
                                if (Double.isNaN(x)) {
                                    System.out.println("m");
                                }
                                dist_new[pId[i][j][a]][pId[i][j][b]] = dist_new[pId[i][j][b]][pId[i][j][a]] = dw[0];
                                weights_new[pId[i][j][a]][pId[i][j][b]] = weights_new[pId[i][j][b]][pId[i][j][a]] = dw[1];
                                sumWeights[pId[i][j][a]] += dw[1];
                                sumWeights[pId[i][j][b]] += dw[1];
                                double sij = 0;
                                for (int k = 0; k < d; k++) {
                                    sij += Math.pow((coord[pId[i][j][a]][k] - coord[pId[i][j][b]][k]), 2);
                                }
                                if (sij != 0) {
                                    sij = dist_new[pId[i][j][a]][pId[i][j][b]] / Math.sqrt(sij);
                                }


                                for (int k = 0; k < d; k++) {
                                    positions[pId[i][j][a]][k] += weights_new[pId[i][j][a]][pId[i][j][b]] * (coord[pId[i][j][b]][k] + sij * (coord[pId[i][j][a]][k] - coord[pId[i][j][b]][k]));
                                    positions[pId[i][j][b]][k] += weights_new[pId[i][j][a]][pId[i][j][b]] * (coord[pId[i][j][a]][k] + sij * (coord[pId[i][j][b]][k] - coord[pId[i][j][a]][k]));
                                    if (Double.isNaN(positions[pId[i][j][a]][k])) {
                                        // positions[pId[i][j][a]][k] = 0.1; 
                                        //System.out.println("m");
                                    }
                                }
                            }
                        }

                    }



                    //pId[i + 1][j][b]
                    if (i + 1 < numBinsX) {
                        for (int a = 0; a < pId[i][j].length; a++) {
                            for (int b = 0; b < pId[i + 1][j].length; b++) {
                                double x = dist(coord[pId[i][j][a]], coord[pId[i + 1][j][b]]);
                                if (!g.isNeighbor(pId[i][j][a], pId[i + 1][j][b])) {
                                    double[] dw = s.parabola(x, false);
                                    if (Double.isNaN(x)) {
                                        System.out.println("m");
                                    }
//                                    if (Double.isNaN(dw[0])) {
//                                        System.out.println("m");
//                                    }
                                    dist_new[pId[i][j][a]][pId[i + 1][j][b]] = dist_new[pId[i + 1][j][b]][pId[i][j][a]] = dw[0];
                                    weights_new[pId[i][j][a]][pId[i + 1][j][b]] = weights_new[pId[i + 1][j][b]][pId[i][j][a]] = dw[1];
                                    sumWeights[pId[i][j][a]] += dw[1];
                                    sumWeights[pId[i + 1][j][b]] += dw[1];
                                    double sij = 0;
                                    for (int k = 0; k < d; k++) {
                                        sij += Math.pow((coord[pId[i][j][a]][k] - coord[pId[i + 1][j][b]][k]), 2);
                                    }
                                    if (sij != 0) {
                                        sij = dist_new[pId[i][j][a]][pId[i + 1][j][b]] / Math.sqrt(sij);
                                    }

                                    for (int k = 0; k < d; k++) {
                                        positions[pId[i][j][a]][k] += weights_new[pId[i][j][a]][pId[i + 1][j][b]] * (coord[pId[i + 1][j][b]][k] + sij * (coord[pId[i][j][a]][k] - coord[pId[i + 1][j][b]][k]));
                                        positions[pId[i + 1][j][b]][k] += weights_new[pId[i][j][a]][pId[i + 1][j][b]] * (coord[pId[i][j][a]][k] + sij * (coord[pId[i + 1][j][b]][k] - coord[pId[i][j][a]][k]));
                                        if (Double.isNaN(positions[pId[i][j][a]][k])) {
                                            // positions[pId[i][j][a]][k] = 0.1; 
                                            //System.out.println("m");
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //pId[i2][j + 1][b]
                    if (j + 1 < numBinsY) {
                        for (int i2 = Math.max(0, i - 1); i2 < Math.min(numBinsX, i + 2); i2++) {
                            for (int a = 0; a < pId[i][j].length; a++) {
                                for (int b = 0; b < pId[i2][j + 1].length; b++) {
                                    double x = dist(coord[pId[i][j][a]], coord[pId[i2][j + 1][b]]);
                                    if (!g.isNeighbor(pId[i][j][a], pId[i2][j + 1][b])) {
                                        double[] dw = s.parabola(x, false);
                                        if (Double.isNaN(x)) {
                                            System.out.println("m");
                                        }
//                                        if (Double.isNaN(dw[0])) {
//                                            System.out.println("m");
//                                        }
                                        dist_new[pId[i][j][a]][pId[i2][j + 1][b]] = dist_new[pId[i2][j + 1][b]][pId[i][j][a]] = dw[0];
                                        weights_new[pId[i][j][a]][pId[i2][j + 1][b]] = weights_new[pId[i2][j + 1][b]][pId[i][j][a]] = dw[1];
                                        sumWeights[pId[i][j][a]] += dw[1];
                                         sumWeights[pId[i2][j + 1][b]] += dw[1];
                                        double sij = 0;
                                        for (int k = 0; k < d; k++) {
                                            sij += Math.pow((coord[pId[i][j][a]][k] - coord[pId[i2][j + 1][b]][k]), 2);
                                        }
                                        if (sij != 0) {
                                            sij = dist_new[pId[i][j][a]][pId[i2][j + 1][b]] / Math.sqrt(sij);
                                        }
                                        for (int k = 0; k < d; k++) {
                                            positions[pId[i][j][a]][k] += weights_new[pId[i][j][a]][pId[i2][j + 1][b]] * (coord[pId[i2][j + 1][b]][k] + sij * (coord[pId[i][j][a]][k] - coord[pId[i2][j + 1][b]][k]));
                                            positions[pId[i2][j + 1][b]][k] += weights_new[pId[i][j][a]][pId[i2][j + 1][b]] * (coord[pId[i][j][a]][k] + sij * (coord[pId[i2][j + 1][b]][k] - coord[pId[i][j][a]][k]));
                                            if (Double.isNaN(positions[pId[i][j][a]][k])) {
                                                //System.out.println("m");
                                                // positions[pId[i][j][a]][k] = 0.1; 
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }



            for (int i = 0; i < g.getVertexCount(); i++) {
                // if(i == 367)
                //  System.out.println("n of 367");
                Collection<Integer> n = g.getNeighbors(i);
                for (int j : n) {
                    if (i < j) {
                        double x = dist(coord[i], coord[j]);
                        double[] dw = s.parabola(x, true);
                        if (Double.isNaN(x)) {
                            System.out.println("m");
                        }
//                        if (Double.isNaN(dw[0])) {
//                            System.out.println("m");
//                        }
                        dist_new[i][j] = dist_new[j][i] = dw[0];
                        weights_new[i][j] = weights_new[j][i] = dw[1];
//                        dwUpdated[i][j]++;
//                        dwUpdated[j][i]++;
                        sumWeights[i] += dw[1];
                        sumWeights[j] += dw[1];
                        double sij = 0;
                        for (int k = 0; k < d; k++) {
                            sij += Math.pow((coord[i][k] - coord[j][k]), 2);
                        }
                        if (sij != 0) {
                            sij = dist_new[i][j] / Math.sqrt(sij);
                        }
                        for (int k = 0; k < d; k++) {
                            positions[i][k] += weights_new[i][j] * (coord[j][k] + sij * (coord[i][k] - coord[j][k]));
                            positions[j][k] += weights_new[i][j] * (coord[i][k] + sij * (coord[j][k] - coord[i][k]));
                            if (Double.isNaN(positions[i][k])) {
                                // positions[i][k] = 0.1;
                                //System.out.println("m");
                            }
                        }

                    }
                }

            }


//            //DEBUG: ok, when all objects in one grid cell
//            boolean allUpdated = true;
//            int count = 0;
//            for (int i = 0; i < n; i++) {
//                for (int j = 0; j < n; j++) {
//                    if (i != j && dwUpdated[i][j] != 1) {
//                        allUpdated = false;
//                        count++;
//                    }
//                }
//            }
//
//            //DEBUG
            //   System.out.println(allUpdated + " " + count);

//            for (int i = 0; i < g.getVertexCount(); i++) {
//                for (int j = 0; j < i; j++) {
//                    double x = dist(coord[i], coord[j]);
//
//
//                    //double aktW = ((2 * s.mu - 2 * x) * (s.a - s.b)) / (Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));
//
//                    double[] dw = s.parabola(x, graph.isNeighbor(i, j));
//                    dist_new[i][j] = dist_new[j][i] = dw[0];
//                    weights_new[i][j] = weights_new[j][i] = dw[1];
//                    if (weights_new[i][j] < minWeight) {
//                        weights_new[i][j] = weights_new[j][i] = minWeight;
//                        correctionCounter++;
//                    }
//                    sumWeights_new += weights_new[i][j];
//                    sumWeights_old += weights[i][j];
//
//                }
//            }
            System.out.println(loop + " " + correctionCounter + " " + m + " " + s.sigma + " " + aktCost);
            //Test
            //  weights_new = scale(weights_new);
            //dist_new = scale(dist_new);
            //DataUtils du = new DataUtils();

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
//            WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
//            wp.iterate();
//            double[][] bla = wp.getPositions();
//


            //double[][] bla = MDSJ.stressMinimization(dist_new, weights_new, 2);

            // double[][] coord_new = new Matrix(bla).transpose().getArrayCopy();
            //gk.setCoord(coord_new);


            // System.out.println("m");

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


////            //DEBUG
////            //check if positions NaN here
//            boolean okBeforeScaling = true;
//            for (int i = 0; i < n; i++) {
//                for (int j = 0; j < d; j++) {
//                    if (Double.isNaN(positions[i][j])) {
//                        okBeforeScaling = false;
//                    }
//                }
//            }
//
//            boolean okAfterScaling = true;
            //copy coord, distance and weights

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    //sumWeights[i] = Math.max(sumWeights[i], 1E-9);
                    if (!(positions[i][j] == 0 && sumWeights[i] == 0)) {
                        //   System.out.println("m");
                        coord[i][j] = positions[i][j] / sumWeights[i];
                    }
                    if (Double.isNaN(coord[i][j])) {
                        //okAfterScaling = false;
                    }

                }
            }
            //        System.out.println(loop + " okBeforeScaling: " + okBeforeScaling + " okAfterScaling: " + okAfterScaling);
//             coord = du.colScaleData(coord);





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
            S -=
                    (Se - z) / (S1e - (S + 2) * (Se - z) / (2 * S + 2));
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

    private double gridThreshold() {
        double xiOld = 1.0;
        double xiNew = 0.0;
        boolean converged = false;
        double convConst = 1E-6;
        double xi = 1.0;
        while (!converged) {
            double xiTerm = Math.pow(tau, 2) / Math.pow(xiOld, 2);
            double product = 4 * Math.PI * xiTerm * Math.pow(Math.log(2), 2) * Math.pow(sigma, 4);
            xiNew = Math.sqrt(-Math.log(product));
            if (Math.abs(xiOld - xiNew) < convConst) {
                converged = true;
            }
            xiOld = xiNew;
        }
        double dist = Math.sqrt(2.0) * sigma * xiNew + mu;
        return dist;
    }
}
