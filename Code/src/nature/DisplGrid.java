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
public class DisplGrid extends AbstractLayout<Integer, Integer> implements IterativeContext {

    private int currentIteration;
    // private int loop;
    private double bestCost; //overall best cost
    private double lastCost = Double.MAX_VALUE; //cost in previous interation
    private int d = 2;
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
    double[][] positions;
    double[][] coord;
    double[][] groundTruth; //ground truth coordinates
    boolean isoInit = false;
    boolean mdsInit = false;
    boolean groundTruthInit = false;
    boolean randomInit = true;
    int convergeCounter = 0;
    boolean varianceFree = true;
    SimpleSigmoid s;
    PairOwn[] p;
    Graph g;
    GraphCompression gk;
    // int d;
    int n;
    int m;
    //Random r;
    //double factor; //how much to keep of old coordinate
    Dimension dd; //size of display
    int numBins;
    static int maxIteration = 500;

    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    public void setXydata(Point2D[] xydata) {
        this.xydata = xydata;
        scaleToDisplaySize();
    }

    public DisplGrid(Graph<Integer, Integer> g) {
        super(g);
        this.g = g;
        // = new Random(20);
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

    public DisplGrid(Graph<Integer, Integer> g, double factor) {
        super(g);
        // r = new Random(20);
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

    //27.4.2015: quadratic space complexity still requried for initialization
    public void initialize() {
        currentIteration = 0;
        // = 0;
        dd = getSize();

        vertices = new int[n];
        for (int i = 0; i < n; i++) {
            vertices[i] = i;
        }

        coord = new double[n][d];
        positions = new double[n][d];
        //init coord with isomap
        if (isoInit) {
            double[][] dist = pathdist();
            coord = new Matrix(MDSJ.classicalScaling(dist)).transpose().getArrayCopy();
        }
        if (groundTruthInit) {
            coord = groundTruth;
            double[][] distances = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
                }
            }

        }

        if (randomInit) {
            Random rr = new Random(10);
            for (int i = 0; i < coord.length; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = rr.nextDouble();
                }
            }
        }

        if (mdsInit) {
            double[][] dd = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    if (g.isNeighbor(i, j)) {
                        dd[i][j] = dd[j][i] = 1.0;
                    } else {
                        dd[i][j] = dd[j][i] = 2.0;
                    }
                }
                coord = new Matrix(MDSJ.stressMinimization(dd)).transpose().getArrayCopy();
            }

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
        gk.setCoord(coord);
        // = 0;
        bestCost = gk.mdlFunctionSimpleSigmoid();
        paramCost = gk.paramCostBic();
        mu = gk.s.mu;
        sigma = gk.s.sigma;
        tau = 1E-6;
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

    private synchronized void update() {
        try {


            DataUtils du = new DataUtils();


            double cutOff = gridThreshold();
            //cutOff = 2.0;
            
            //cutOff = 100;
            double[][] minMax = du.minMax(coord);
            int numBinsX = (int) Math.ceil((minMax[0][1] - minMax[0][0]) / cutOff);
            int numBinsY = (int) Math.ceil((minMax[1][1] - minMax[1][0]) / cutOff);
            // System.out.println(cutOff + " " + numBinsX + " " + numBinsY);


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

if(currentIteration % 30 == 0){
            //compute sigmoid
            p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), g.isNeighbor(i, j));
                    //DEBUG
                    if (Double.isNaN(p[i * (i - 1) / 2 + j].dist)) {
                        //    System.out.println("dist NaN");
                    }
                    //DEBUG
                }
            }

             s = new SimpleSigmoid(p);
            mu = s.mu;
            sigma = s.sigma;



}
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

                                double[] dw = s.parabola(x, false);
                                if (Double.isNaN(x)) {
                                    System.out.println("m");
                                }
                                // dist_new[pId[i][j][a]][pId[i][j][b]] = dist_new[pId[i][j][b]][pId[i][j][a]] = dw[0];
                                // weights_new[pId[i][j][a]][pId[i][j][b]] = weights_new[pId[i][j][b]][pId[i][j][a]] = dw[1];
                                sumWeights[pId[i][j][a]] += dw[1];
                                sumWeights[pId[i][j][b]] += dw[1];
                                double sij = 0;
                                for (int k = 0; k < d; k++) {
                                    sij += Math.pow((coord[pId[i][j][a]][k] - coord[pId[i][j][b]][k]), 2);
                                }
                                if (sij != 0) {
                                    sij = dw[0] / Math.sqrt(sij);
                                }


                                for (int k = 0; k < d; k++) {
                                    positions[pId[i][j][a]][k] += dw[1] * (coord[pId[i][j][b]][k] + sij * (coord[pId[i][j][a]][k] - coord[pId[i][j][b]][k]));
                                    positions[pId[i][j][b]][k] += dw[1] * (coord[pId[i][j][a]][k] + sij * (coord[pId[i][j][b]][k] - coord[pId[i][j][a]][k]));
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
                                    // dist_new[pId[i][j][a]][pId[i + 1][j][b]] = dist_new[pId[i + 1][j][b]][pId[i][j][a]] = dw[0];
                                    // weights_new[pId[i][j][a]][pId[i + 1][j][b]] = weights_new[pId[i + 1][j][b]][pId[i][j][a]] = dw[1];
                                    sumWeights[pId[i][j][a]] += dw[1];
                                    sumWeights[pId[i + 1][j][b]] += dw[1];
                                    double sij = 0;
                                    for (int k = 0; k < d; k++) {
                                        sij += Math.pow((coord[pId[i][j][a]][k] - coord[pId[i + 1][j][b]][k]), 2);
                                    }
                                    if (sij != 0) {
                                        sij = dw[0] / Math.sqrt(sij);
                                    }

                                    for (int k = 0; k < d; k++) {
                                        positions[pId[i][j][a]][k] += dw[1] * (coord[pId[i + 1][j][b]][k] + sij * (coord[pId[i][j][a]][k] - coord[pId[i + 1][j][b]][k]));
                                        positions[pId[i + 1][j][b]][k] += dw[1] * (coord[pId[i][j][a]][k] + sij * (coord[pId[i + 1][j][b]][k] - coord[pId[i][j][a]][k]));
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
                                        // dist_new[pId[i][j][a]][pId[i2][j + 1][b]] = dist_new[pId[i2][j + 1][b]][pId[i][j][a]] = dw[0];
                                        // weights_new[pId[i][j][a]][pId[i2][j + 1][b]] = weights_new[pId[i2][j + 1][b]][pId[i][j][a]] = dw[1];
                                        sumWeights[pId[i][j][a]] += dw[1];
                                        sumWeights[pId[i2][j + 1][b]] += dw[1];
                                        double sij = 0;
                                        for (int k = 0; k < d; k++) {
                                            sij += Math.pow((coord[pId[i][j][a]][k] - coord[pId[i2][j + 1][b]][k]), 2);
                                        }
                                        if (sij != 0) {
                                            sij = dw[0] / Math.sqrt(sij);
                                        }
                                        for (int k = 0; k < d; k++) {
                                            positions[pId[i][j][a]][k] += dw[1] * (coord[pId[i2][j + 1][b]][k] + sij * (coord[pId[i][j][a]][k] - coord[pId[i2][j + 1][b]][k]));
                                            positions[pId[i2][j + 1][b]][k] += dw[1] * (coord[pId[i][j][a]][k] + sij * (coord[pId[i2][j + 1][b]][k] - coord[pId[i][j][a]][k]));
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

                        sumWeights[i] += dw[1];
                        sumWeights[j] += dw[1];
                        double sij = 0;
                        for (int k = 0; k < d; k++) {
                            sij += Math.pow((coord[i][k] - coord[j][k]), 2);
                        }
                        if (sij != 0) {
                            sij = dw[0] / Math.sqrt(sij);
                        }
                        for (int k = 0; k < d; k++) {
                            positions[i][k] += dw[1] * (coord[j][k] + sij * (coord[i][k] - coord[j][k]));
                            positions[j][k] += dw[1] * (coord[i][k] + sij * (coord[j][k] - coord[i][k]));
                            if (Double.isNaN(positions[i][k])) {
                                // positions[i][k] = 0.1;
                                //System.out.println("m");
                            }
                        }

                    }
                }

            }



            System.out.println(currentIteration + " " + correctionCounter + " " + m + " " + s.sigma + " " + aktCost 
            + " " + cutOff + " " + numBinsX + " " + numBinsY);
            //Test
            //  weights_new = scale(weights_new);
            //dist_new = scale(dist_new);
            //DataUtils du = new DataUtils();

            String name = "results_" + currentIteration;


            String name1 = "dist_" + currentIteration;
            String name2 = "weights_" + currentIteration;
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


            String name3 = "coord_" + currentIteration;
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
//                bestloop = loop;
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

            // loop++;
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
