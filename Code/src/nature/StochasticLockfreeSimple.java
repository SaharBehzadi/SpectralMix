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
public class StochasticLockfreeSimple extends AbstractLayout<Integer, Integer> implements IterativeContext {

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
    double[][] coordOld; //relevant
    double[][] coord; //relevant
    double[] nodeWeights; //relevant
    double[][] weights;
    double[] sumWeights;
    EdgeUpdate[] e;
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
    // int d;
    int n;
    int m;
    //double factor; //how much to keep of old coordinate
    Dimension dd; //size of display
    private static double PRECISION = 1E-12;
    static int maxIteration = 500;

    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    public void setXydata(Point2D[] xydata) {
        this.xydata = xydata;
        scaleToDisplaySize();
    }

    public StochasticLockfreeSimple(Graph<Integer, Integer> g) {
        super(g);
        this.g = g;
        r = new Random(1);
        d = 2;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        coord = new double[n][d];
        coordOld = new double[n][d];
        e = new EdgeUpdate[g.getEdgeCount()];
        nodeWeights = new double[n];

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

    public StochasticLockfreeSimple(Graph<Integer, Integer> g, double factor) {
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
        currentIteration = 0;
        cheatCounter = 0;
        cheatVar = false;
        vertices = new int[n];
        for (int i = 0; i < n; i++) {
            vertices[i] = i;
        }

        coord = new double[n][d];
        coordOld = new double[n][d];
        sumWeights = new double[n];
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
            for (int i = 0; i < coordOld.length; i++) {
                for (int j = 0; j < d; j++) {
                    coordOld[i][j] = r.nextDouble();
                }
            }
            //distances = new double[n][n];
//            for (int i = 0; i < n; i++) {
//                for (int j = 0; j < i; j++) {
//                    distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
//                }
//            }
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
        e = new EdgeUpdate[g.getEdgeCount()];
        computeSigmoid();
        s.sigma = Math.min(s.sigma, 1.0);
        //one round of sampling-based update; initializes e and coord
        double[][] dd = new double[n][n]; //for more memory-efficient implementation this needs to be changed
        double[][] ww = new double[n][n];
        for (int i = 0; i < g.getEdgeCount(); i++) {
            Pair<Integer> akt = g.getEndpoints(i);

            double x = dist(coordOld[akt.getFirst()], coordOld[akt.getSecond()]);
            double[] dw = s.parabola(x, true);
            nodeWeights[akt.getFirst()] += dw[1];
            nodeWeights[akt.getSecond()] += dw[1];
            dd[akt.getFirst()][akt.getSecond()] = dd[akt.getSecond()][akt.getFirst()] = dw[0];
            ww[akt.getFirst()][akt.getSecond()] = ww[akt.getSecond()][akt.getFirst()] = dw[1];
            int notEdge1 = r.nextInt(n);
            while (g.isNeighbor(akt.getFirst(), notEdge1) || akt.getFirst() == notEdge1) {
                notEdge1 = r.nextInt(n);
            }
            x = dist(coordOld[akt.getFirst()], coordOld[notEdge1]);
            dw = s.parabola(x, false);
            dd[akt.getFirst()][notEdge1] = dd[notEdge1][akt.getFirst()] = dw[0];
            ww[akt.getFirst()][notEdge1] = ww[notEdge1][akt.getFirst()] = dw[1];
            nodeWeights[akt.getFirst()] += dw[1];
            nodeWeights[notEdge1] += dw[1];
            int notEdge2 = r.nextInt(n);

            while (g.isNeighbor(akt.getSecond(), notEdge2) || akt.getSecond() == notEdge2) {
                notEdge2 = r.nextInt(n);
            }
            x = dist(coordOld[akt.getSecond()], coordOld[notEdge2]);
            dw = s.parabola(x, false);
            dd[akt.getSecond()][notEdge2] = dd[notEdge2][akt.getSecond()] = dw[0];
            ww[akt.getSecond()][notEdge2] = ww[notEdge2][akt.getSecond()] = dw[1];
            nodeWeights[akt.getSecond()] += dw[1];
            nodeWeights[notEdge2] += dw[1];
            e[i] = new EdgeUpdate(akt.getFirst(), akt.getSecond(), notEdge1, notEdge2);

        }

        WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dd, new Matrix(coordOld).transpose().getArrayCopy(), ww);
        wp.iterate();
        double[][] bla = wp.getPositions();

        coord = new Matrix(bla).transpose().getArrayCopy();

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
            // System.out.println(currentIteration + " " + s.sigma);
            }
            s.sigma = Math.min(s.sigma, 1.0);

            double aktCost = s.costAllPairs(p);

            System.out.println(loop + " " + s.mu + " " + s.sigma + " " + aktCost);
            // if(cheatVar)
            //  s.sigma = Math.max(s.sigma, 5.0);
            // System.out.println("variance before update: " + s.sigma);
            //compute new weights and dists by setting derivatives of MDS error function and cost function equal
            int correctionCounter = 0;
            r = new Random(currentIteration+1);
            for (int i = 0; i < g.getEdgeCount(); i++) {

                Pair<Integer> akt = g.getEndpoints(i);
                updateE(akt, i);

            }

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    // bestDb[i][j] = coord[i][j];
                }
                xydata[i].setLocation(coord[i][0], coord[i][1]);
            }
            scaleToDisplaySize();
            currentIteration++;
        } catch (ConcurrentModificationException cme) {
        }

    }

//    private void updateE(Pair<Integer> edge, boolean connected) {
//        double x = dist(coord[edge.getFirst()], coord[edge.getSecond()]);
//        double[] dw = s.parabola(x, connected);
//        double[] coordNewFirst = new double[d];
//        double[] coordNewSecond = new double[d];
//        double nodeWeightsOld
//        
//        
//        double nodeWeightNewFirst = nodeWeights[edge.getFirst()] + dw[1] - edgeWeights[index_ij];
//        double nodeWeightNewSecond = nodeWeights[edge.getSecond()] + dw[1] - edgeWeights[index_ij];
//        double sij = 0;
//        for (int k = 0; k < d; k++) {
//            sij += Math.pow((coord[edge.getFirst()][k] - coord[edge.getSecond()][k]), 2);
//        }
//        if (sij != 0) {
//            sij = dw[0] / Math.sqrt(sij);
//        }
//        double[] zijNew = new double[d];
//        double[] zjiNew = new double[d];
//        for (int k = 0; k < d; k++) {
//            zijNew[k] = dw[1] * (coord[edge.getSecond()][k] + sij * (coord[edge.getFirst()][k] - coord[edge.getSecond()][k]));
//            coordNewFirst[k] = (nodeWeights[edge.getFirst()] * coord[edge.getFirst()][k] - z[edge.getFirst()][edge.getSecond()][k] + zijNew[k]) / nodeWeightNewFirst;
//        }
//        for (int k = 0; k < d; k++) {
//            zjiNew[k] = dw[1] * (coord[edge.getFirst()][k] + sij * (coord[edge.getSecond()][k] - coord[edge.getFirst()][k]));
//            coordNewSecond[k] = (nodeWeights[edge.getSecond()] * coord[edge.getSecond()][k] - z[edge.getSecond()][edge.getFirst()][k] + zjiNew[k]) / nodeWeightNewSecond;
//        }
//        edgeWeights[index_ij] = dw[1];
//        nodeWeights[edge.getFirst()] = nodeWeightNewFirst;
//        nodeWeights[edge.getSecond()] = nodeWeightNewSecond;
//        coord[edge.getFirst()] = coordNewFirst;
//        coord[edge.getSecond()] = coordNewSecond;
//        z[edge.getFirst()][edge.getSecond()] = zijNew;
//        z[edge.getSecond()][edge.getFirst()] = zjiNew;
//
//    }
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

    //changes: coord and coordOld of both endpoints, nodeWeights of both endpoints
    private void updateE(Pair<Integer> akt, int edgeIndex) {
        //sample new pair of not-edges
        int notEdge1 = r.nextInt(n);
        while (g.isNeighbor(akt.getFirst(), notEdge1)) {
            notEdge1 = r.nextInt(n);
        }
        int notEdge2 = r.nextInt(n);
        while (g.isNeighbor(akt.getSecond(), notEdge2)) {
            notEdge2 = r.nextInt(n);
        }
        //Variables to be modified: 6 Points, their coords and their weights are modified by one edge
        //the two endpoints
        double[] coordNewFirst = new double[d];
        double[] coordNewSecond = new double[d];
        double weightNewFirst = nodeWeights[akt.getFirst()];
        double weightNewSecond = nodeWeights[akt.getSecond()];
        //the four not edges
        double[] coordNewNotEdge1Old = new double[d];
        double[] coordNewNotEdge2Old = new double[d];
        double[] coordNewNotEdge1New = new double[d];
        double[] coordNewNotEdge2New = new double[d];
        double weightNewNotEdge1Old = nodeWeights[e[edgeIndex].pNotEdgeI];
        double weightNewNotEdge2Old = nodeWeights[e[edgeIndex].pNotEdgeJ];
        double weightNewNotEdge1New = nodeWeights[notEdge1];
        double weightNewNotEdge2New = nodeWeights[notEdge2];

        //update the weights and coordinates
        //the mutual influence of the two endpoints       
        double dd = dist(coord[akt.getFirst()], coord[akt.getSecond()]);
        double[] dw = s.parabola(dd, true);
        double dOld = dist(coordOld[akt.getFirst()], coordOld[akt.getSecond()]);
        double[] dwOld = s.parabola(dOld, true);
        double[] zOld = zij(dwOld[0], dwOld[1], akt.getFirst(), akt.getSecond(), true);
        double[] zNew = zij(dw[0], dw[1], akt.getFirst(), akt.getSecond(), false);
        weightNewFirst -= dwOld[1];
        weightNewSecond -= dwOld[1];
        weightNewFirst += dw[1];
        weightNewSecond += dw[1];
        for (int k = 0; k < d; k++) {
            coordNewFirst[k] = (nodeWeights[akt.getFirst()] * coord[akt.getFirst()][k] - zOld[k] + zNew[k]) / weightNewFirst;
        }
        zOld = zij(dwOld[0], dwOld[1], akt.getSecond(), akt.getFirst(), true);
        zNew = zij(dw[0], dw[1], akt.getSecond(), akt.getFirst(), false);
        for (int k = 0; k < d; k++) {
            coordNewSecond[k] = (nodeWeights[akt.getSecond()] * coord[akt.getSecond()][k] - zOld[k] + zNew[k]) / weightNewSecond;
        }
        //the mutual influence of the two not-edges: Probably: consider only the influence of the not edges on the corresponding endpoint and not of the endpoint on the not edge
        dd = dist(coord[akt.getFirst()], coord[notEdge1]);
        dw = s.parabola(dd, false);
        dOld = dist(coordOld[akt.getFirst()], coordOld[e[edgeIndex].pNotEdgeI]);
        dwOld = s.parabola(dOld, false);
        zOld = zij(dwOld[0], dwOld[1], akt.getFirst(), e[edgeIndex].pNotEdgeI, true);
        zNew = zij(dwOld[0], dwOld[1], akt.getFirst(), notEdge1, false);
        weightNewFirst -= dwOld[1];
        weightNewFirst += dw[1];
        weightNewNotEdge1Old -= dwOld[1]; //make sure that nodeWeights > 0
        weightNewNotEdge1Old = Math.max(weightNewNotEdge1Old, 0.0);
        weightNewNotEdge1New += dw[1];
        for (int k = 0; k < d; k++) {
            coordNewFirst[k] = (nodeWeights[akt.getFirst()] * coord[akt.getFirst()][k] - zOld[k] + zNew[k]) / weightNewFirst;
        }
        zOld = zij(dwOld[0], dwOld[1], e[edgeIndex].pNotEdgeI, akt.getFirst(), true);
        zNew = zij(dw[0], dw[1], notEdge1, akt.getFirst(), false);
        for (int k = 0; k < d; k++) {
            coordNewNotEdge1Old[k] = (nodeWeights[e[edgeIndex].pNotEdgeI] * coord[e[edgeIndex].pNotEdgeI][k] - zOld[k]) / weightNewNotEdge1Old;
        }
        for (int k = 0; k < d; k++) {
            coordNewNotEdge1New[k] = (nodeWeights[notEdge1] * coord[notEdge1][k] + zNew[k]) / weightNewNotEdge1New;
        }
        
        dd = dist(coord[akt.getSecond()], coord[notEdge2]);
        dw = s.parabola(dd, false);
        dOld = dist(coordOld[akt.getSecond()], coordOld[e[edgeIndex].pNotEdgeJ]);
        dwOld = s.parabola(dOld, false);
        zOld = zij(dwOld[0], dwOld[1], akt.getSecond(), e[edgeIndex].pNotEdgeJ, true);
        zNew = zij(dwOld[0], dwOld[1], akt.getSecond(), notEdge2, false);
        weightNewSecond -= dwOld[1];
        weightNewSecond += dw[1];
        weightNewNotEdge2Old -= dwOld[1]; //make sure that nodeWeights > 0
        weightNewNotEdge2Old = Math.max(weightNewNotEdge2Old, 0.0);
        weightNewNotEdge2New += dw[1];
        for (int k = 0; k < d; k++) {
            coordNewSecond[k] = (nodeWeights[akt.getSecond()] * coord[akt.getSecond()][k] - zOld[k] + zNew[k]) / weightNewSecond;
        }
        zOld = zij(dwOld[0], dwOld[1], e[edgeIndex].pNotEdgeJ, akt.getSecond(), true);
        zNew = zij(dw[0], dw[1], notEdge2, akt.getSecond(), false);
        for (int k = 0; k < d; k++) {
            coordNewNotEdge2Old[k] = (nodeWeights[e[edgeIndex].pNotEdgeJ] * coord[e[edgeIndex].pNotEdgeJ][k] - zOld[k]) / weightNewNotEdge2Old;
        }
        for (int k = 0; k < d; k++) {
            coordNewNotEdge1New[k] = (nodeWeights[notEdge1] * coord[notEdge1][k] + zNew[k]) / weightNewNotEdge1New;
        }
        
        
        //write the updates
        e[edgeIndex].pNotEdgeI = notEdge1;
        e[edgeIndex].pNotEdgeJ = notEdge2;
        writeUpdates(akt.getFirst(), coordNewFirst, weightNewFirst);
        writeUpdates(akt.getSecond(), coordNewSecond, weightNewSecond);
        writeUpdates(e[edgeIndex].endpointI, coordNewNotEdge1Old, weightNewNotEdge1Old);
        writeUpdates(e[edgeIndex].endpointJ, coordNewNotEdge2Old, weightNewNotEdge2Old);
        writeUpdates(notEdge1, coordNewNotEdge1New, weightNewNotEdge1New);
        writeUpdates(notEdge2, coordNewNotEdge2New, weightNewNotEdge2New);

    }

    private void writeUpdates(int index, double[] coordNew, double weightNew) {
        nodeWeights[index] = weightNew;
        double[] coordOldR = new double[d];
        for (int i = 0; i < d; i++) {
            coordOldR[i] = coord[index][i];
        }
        coordOld[index] = coordOldR;
        coord[index] = coordNew;
        nodeWeights[index] = weightNew;

    }

    private double[] zij(double dij, double wij, int i, int j, boolean old) {
        double[] res = new double[d];
        double sij = 0.0;
        for (int k = 0; k < d; k++) {
            if (old) {
                sij += Math.pow((coordOld[i][k] - coordOld[j][k]), 2);
            } else {
                sij += Math.pow((coord[i][k] - coord[j][k]), 2);
            }
        }
        if (sij != 0) {
            sij = dij / Math.sqrt(sij);
        }
        for (int k = 0; k < d; k++) {
            if (old) {
                res[k] = wij * (coordOld[j][k] + sij * (coordOld[i][k] - coordOld[j][k]));
            } else {
                res[k] = wij * (coord[j][k] + sij * (coord[i][k] - coord[j][k]));
            }
        }
        return res;
    }
}
