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
public class IncrementalOwnDispl extends AbstractLayout<Integer, Integer> implements IterativeContext {

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
    boolean isoInit = false;
    boolean mdsInit = false;
    boolean groundTruthInit = false;
    boolean randomInit = true;
    int convergeCounter = 0;
    boolean varianceFree = false;
    Graph g;
    GraphCompression gk;
    SimpleSigmoid s;
    Random r;
    double[] nodeWeights;
    double[] edgeWeights;
    double[][][] z; //n x n x d
    // int d;
    int n;
    int m;
    //double factor; //how much to keep of old coordinate
    Dimension dd; //size of display
    private static double PRECISION = 1E-12;
    static int maxIteration = 5000;
    boolean verbose = true;

    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    public void setXydata(Point2D[] xydata) {
        this.xydata = xydata;
        scaleToDisplaySize();
    }

    public IncrementalOwnDispl(Graph<Integer, Integer> g) {
        super(g);
        this.g = g;
        r = new Random(1);
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

    public IncrementalOwnDispl(Graph<Integer, Integer> g, double factor) {
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
        r = new Random(1);
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
            //init with one full iteration of majorization
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        for (int k = 0; k < d; k++) {
                            coord[i][k] = r.nextDouble();
                        }
                    }
                }
            }
            computeSigmoid();
            s.sigma = Math.min(s.sigma, 1.0);
            nodeWeights = new double[n];
            edgeWeights = new double[m];
            z = new double[n][n][d];
            double[][] weights_new = new double[n][n];
            double[][] dist_new = new double[n][n];
            DataUtils du = new DataUtils();

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    double[] dw = s.parabola(dist(coord[i], coord[j]), g.isNeighbor(i, j));
                    weights_new[i][j] = weights_new[j][i] = dw[1];
                    dist_new[i][j] = dist_new[j][i] = dw[0];
                    int index_ij = du.getIndex(i, j, n);
                    edgeWeights[index_ij] = dw[1];

                }
            }
            //init zij and nodeweights
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j) {
//                    if(j == 23)
//                        System.out.println("m");
                        nodeWeights[i] += weights_new[i][j];
                        double sij = 0;
                        for (int k = 0; k < d; k++) {
                            sij += Math.pow((coord[i][k] - coord[j][k]), 2);
                        }
                        //   sij = Math.max(sij, 1e-200);
                        if (sij != 0) {
                            sij = dist_new[i][j] / Math.sqrt(sij);
                        }

                        for (int k = 0; k < d; k++) {
                            z[i][j][k] += weights_new[i][j] * (coord[j][k] + sij * (coord[i][k] - coord[j][k]));
//                        if (Double.isNaN(positions[i][k])) {
//                           System.out.println("i " + i + " k: " + k + " j: " + j);
//                        }
                        }
                    }
                }//j
            }

            WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
            wp.iterate();
            double[][] bla = wp.getPositions();

            double[][] coord_new = new Matrix(bla).transpose().getArrayCopy();
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = coord_new[i][j];
                }
            }
//TEST check if coord[0][0] is sum of z[0][:][0]/nodeweights[0]: ok
//            double testSum = 0;
//            for (int i = 0; i < n; i++) {
//                testSum += z[0][i][0];
//            }
//            testSum /= nodeWeights[0];
//            System.out.println("testSum " + testSum + " " + coord[0][0]);

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
//        computeSigmoid();
//        gk.setCoord(coord);
//        loop = 0;
//        bestCost = gk.mdlFunctionSimpleSigmoid();
//        paramCost = gk.paramCostBic();
//
//        System.out.println("pCostsBic: " + paramCost);
//        savedBits = costWithoutEmbedding - (bestCost + paramCost);
//        System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
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

    private synchronized void update() {
        try {
            double[][] weights_new = new double[n][n];
            double sumWeights_new = 0.0;
            double sumWeights_old = 0.0;
            double[][] dist_new = new double[n][n];
            double minWeight = 1E-9;
            double minDist = 1E-6;

//
//            //compute sigmoid
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                }
            }
//            SimpleSigmoid s = new SimpleSigmoid(p);
//            if (s.sigma > 1) {
//                s.sigma = 1;

             if (currentIteration % 100 == 0) {
            computeSigmoid();
            //System.out.println(s.sigma);
             }
            s.sigma = Math.min(s.sigma, 1.0);

            double aktCost = s.costAllPairs(p);

            System.out.println(loop + " " + s.sigma + " Cost=" + aktCost);
            // if(cheatVar)
            //  s.sigma = Math.max(s.sigma, 5.0);
            // System.out.println("variance before update: " + s.sigma);
            //compute new weights and dists by setting derivatives of MDS error function and cost function equal
            int correctionCounter = 0;

            for (int i = 0; i < g.getEdgeCount(); i++) {
                Pair<Integer> akt = g.getEndpoints(i);
                updateE(akt, true);
                int notEdge1 = r.nextInt(n);

                while (g.isNeighbor(akt.getFirst(), notEdge1) || akt.getFirst() == notEdge1) {
                    notEdge1 = r.nextInt(n);
                }
                Pair<Integer> e = new Pair<Integer>(akt.getFirst(), notEdge1);
                updateE(e, false);
                int notEdge2 = r.nextInt(n);

                while (g.isNeighbor(akt.getSecond(), notEdge2) || akt.getSecond() == notEdge2) {
                    notEdge2 = r.nextInt(n);
                }
                e = new Pair<Integer>(akt.getSecond(), notEdge2);
                updateE(e, false);

            }
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
//            WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
//            wp.iterate();
//            double[][] bla = wp.getPositions();
//
//            //double[][] bla = MDSJ.stressMinimization(dist_new, weights_new, 2);
//            double[][] coord_new = new Matrix(bla).transpose().getArrayCopy();
//            gk.setCoord(coord_new);
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
//            for (int i = 0; i < n; i++) {
//                for (int j = 0; j < d; j++) {
//                    coord[i][j] = coord_new[i][j];
//                }
//            }
//
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
        } catch (ConcurrentModificationException cme) {
        }

    }

    public void updateE(Pair<Integer> edge, boolean connected) {
        DataUtils du = new DataUtils();
        double x = dist(coord[edge.getFirst()], coord[edge.getSecond()]);
        double[] dw = s.parabola(x, connected);
        double[] coordNewFirst = new double[d];
        double[] coordNewSecond = new double[d];
        int index_ij = du.getIndex(edge.getFirst(), edge.getSecond(), n);
        double nodeWeightNewFirst = nodeWeights[edge.getFirst()] + dw[1] - edgeWeights[index_ij];
        double nodeWeightNewSecond = nodeWeights[edge.getSecond()] + dw[1] - edgeWeights[index_ij];
        double sij = 0;
        for (int k = 0; k < d; k++) {
            sij += Math.pow((coord[edge.getFirst()][k] - coord[edge.getSecond()][k]), 2);
        }
        if (sij != 0) {
            sij = dw[0] / Math.sqrt(sij);
        }
        double[] zijNew = new double[d];
        double[] zjiNew = new double[d];
        for (int k = 0; k < d; k++) {
            zijNew[k] = dw[1] * (coord[edge.getSecond()][k] + sij * (coord[edge.getFirst()][k] - coord[edge.getSecond()][k]));
            coordNewFirst[k] = (nodeWeights[edge.getFirst()] * coord[edge.getFirst()][k] - z[edge.getFirst()][edge.getSecond()][k] + zijNew[k]) / nodeWeightNewFirst;
        }
        for (int k = 0; k < d; k++) {
            zjiNew[k] = dw[1] * (coord[edge.getFirst()][k] + sij * (coord[edge.getSecond()][k] - coord[edge.getFirst()][k]));
            coordNewSecond[k] = (nodeWeights[edge.getSecond()] * coord[edge.getSecond()][k] - z[edge.getSecond()][edge.getFirst()][k] + zjiNew[k]) / nodeWeightNewSecond;
        }
        edgeWeights[index_ij] = dw[1];
        nodeWeights[edge.getFirst()] = nodeWeightNewFirst;
        nodeWeights[edge.getSecond()] = nodeWeightNewSecond;
        coord[edge.getFirst()] = coordNewFirst;
        coord[edge.getSecond()] = coordNewSecond;
        z[edge.getFirst()][edge.getSecond()] = zijNew;
        z[edge.getSecond()][edge.getFirst()] = zjiNew;
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
