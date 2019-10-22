package nature;

//package edu.uci.ics.jung.algorithms.layout;
import Jama.Matrix;
import edu.uci.ics.jung.algorithms.layout3d.AbstractLayout;
import edu.uci.ics.jung.algorithms.layout.GraphElementAccessor;
import edu.uci.ics.jung.algorithms.layout.util.RandomLocationTransformer;
import edu.uci.ics.jung.algorithms.layout3d.SpringLayout.LengthFunction;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.algorithms.shortestpath.Distance;
import edu.uci.ics.jung.algorithms.util.IterativeContext;
import org.apache.commons.math3.special.Erf;


import java.awt.Dimension;
//import java.awt.geom.Point2D;
import javax.vecmath.Point3f;
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
public class WeightedMajorizationDispl3d extends AbstractLayout<Integer, Integer> implements IterativeContext {

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
    private int d = 3;
    private String status = "OwnLayout";
    private boolean adjustForGravity = true;
    private int[] vertices;
    private Point3f[] xydata;
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
    // int d;
    int n;
    int m;
    Random r;
    //double factor; //how much to keep of old coordinate
    //Dimension dd; //size of display
    private static double PRECISION = 1E-12;
    static int maxIteration = 2000;

    public WeightedMajorizationDispl3d(Graph<Integer, Integer> g) {
        super(g);

        //this(g, UNITLENGTHFUNCTION);
        this.g = g;

        r = new Random(20);
        //this.factor = factor;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;

    }

    public WeightedMajorizationDispl3d(Graph<Integer, Integer> g, double[][] coordGroundTruth) {
        super(g);

        //this(g, UNITLENGTHFUNCTION);
        this.g = g;
        this.groundTruth = coordGroundTruth;
        groundTruthInit = true;
        r = new Random(20);
        //this.factor = factor;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;

    }

    public static final class UnitLengthFunction<E> implements LengthFunction<E> {

        int length;

        public UnitLengthFunction(int length) {



            this.length = length;



        }

        public double getLength(E e) {



            return length;



        }
    }
    public static final LengthFunction UNITLENGTHFUNCTION = new UnitLengthFunction(
            30);

    public double getBestllh() {
        return bestCost;
    }

    public double[][] getCoordinates() {

        return bestDb;
    }

    public void setGroundTruth(double[][] groundTruth) {
        this.groundTruth = groundTruth;
    }

    public WeightedMajorizationDispl3d(Graph<Integer, Integer> g, double factor) {
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
//        dd = getSize();
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
            coord = new Matrix(MDSJ.classicalScaling(distances, 3)).transpose().getArrayCopy();
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
            Random r = new Random();
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

        xydata = new Point3f[n];
        bestDb = new double[n][d];
        for (int i = 0; i < n; i++) {
            xydata[i] = transform(vertices[i]);
            float[] p = new float[3];
            p[0] = (float) coord[i][0];
            p[1] = (float) coord[i][1];
            p[2] = (float) coord[i][2];
            xydata[i] = new Point3f(p);
            bestDb[i][0] = coord[i][0];
            bestDb[i][1] = coord[i][1];
            bestDb[i][1] = coord[i][2];
        }
        gk.setCoord(coord);
        loop = 0;
        bestCost = gk.mdlFunctionSimpleSigmoid();
        paramCost = gk.paramCostBic();

        System.out.println("pCostsBic: " + paramCost);
        savedBits = costWithoutEmbedding - (bestCost + paramCost);
        System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
//        scaleToDisplaySize();

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

    private synchronized void update() {
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
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), g.isNeighbor(i, j));
                }
            }
            SimpleSigmoid s = new SimpleSigmoid(p);
            if (s.sigma > 1) {
                s.sigma = 1;
            }
//            double t = getDistThreshold(20);
//            System.out.println("20 " + t);
//            t = getDistThreshold(10);
//            System.out.println("10 " + t);
//            t = getDistThreshold(5);
//            System.out.println("5 " + t);
//            t = getDistThreshold(3);
//            System.out.println("3 " + t);
//            t = getDistThreshold(1);
//            System.out.println("1 " + t);
//             t = getDistThreshold(0.5);
//            System.out.println("0.5 " + t);
//             t = getDistThreshold(0.3);
//            System.out.println("0.3 " + t);
//            t = getDistThreshold(0.1);
//            System.out.println("0.1 " + t);
//            t = getDistThreshold(0.01);
//            System.out.println("0.01 " + t);
//            System.out.println("m");

//            if (convergeCounter == 0) {
//                s.sigma = 20;
//            }
//
//            if (convergeCounter == 1) {
//                s.sigma = 10;
//                lastCost = Double.MAX_VALUE;
//
//            }
//
//            if (convergeCounter == 2) {
//                s.sigma = 5;
//                lastCost = Double.MAX_VALUE;
//
//            }
//
//            if (convergeCounter == 3) {
//                s.sigma = 1;
//                lastCost = Double.MAX_VALUE;
//            }

//            if (loop < 1000 && s.sigma < 5.0) {
//                s.sigma = 5.0;
//            }
////                if(loop == 22)
////                    System.out.println("m");
//                s = new SimpleSigmoid(1.0);
//            }
//            if (loop > 50) {
//                s = new SimpleSigmoid(p);
//            }

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

                    double[] dw = s.parabola(x, g.isNeighbor(i, j));
                    dist_new[i][j] = dist_new[j][i] = dw[0];
                    weights_new[i][j] = weights_new[j][i] = dw[1];
                    if (weights_new[i][j] < minWeight) {
                        weights_new[i][j] = weights_new[j][i] = minWeight;
                        dist_new[i][j] = dist_new[j][i] = distances[i][j];
                        correctionCounter++;
                    }
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
            if(loop == 2000)
            // du.saveAsMatlab3(dist_new, weights_new, coord_new, "dist", "weight", "coord", name+".mat");
             ea.writeDoubleToMatlab(coord, name3);


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
                float[] pc = new float[3];
                pc[0] = (float) coord[i][0];
                pc[1] = (float) coord[i][1];
                pc[2] = (float) coord[i][2];
                xydata[i] = new Point3f(pc);
            }
//            scaleToDisplaySize();


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
            //();
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

    public void step() {
        //randomShuffle();

        update();
//        scaleToDisplaySize();
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
//    public void scaleToDisplaySize() {
//
//
//        double borderWidth = getSize().getRadius() / 50.0;
//        double minDisplay = -getSize().getRadius() + borderWidth;
//        double maxDisplay = -minDisplay;
//        double displaySpread = maxDisplay - minDisplay;
//
//
//        //scale coordinates between 0 and 1
//        double max_x = -java.lang.Double.MAX_VALUE;
//        double max_y = -java.lang.Double.MAX_VALUE;
//        double max_z = -java.lang.Double.MAX_VALUE;
//        double min_x = -max_x;
//        double min_y = -max_y;
//        double min_z = -max_z;
//
//        for (int i = 0; i < n; i++) {
//
//            if (xydata[i].getX() > max_x) {
//                max_x = xydata[i].getX();
//            }
//            if (xydata[i].getX() < min_x) {
//                min_x = xydata[i].getX();
//            }
//            if (xydata[i].getY() > max_y) {
//                max_y = xydata[i].getY();
//            }
//            if (xydata[i].getY() < min_y) {
//                min_y = xydata[i].getY();
//            }
//
//            if (xydata[i].getY() > max_z) {
//                max_z = xydata[i].getZ();
//            }
//            if (xydata[i].getY() < min_z) {
//                min_z = xydata[i].getZ();
//            }
//
//        }
//
//        for (int i = 0; i < xydata.length; i++) {
//            double x_new = minDisplay + (((xydata[i].getX() - min_x) / (max_x - min_x)) * 0.99 * displaySpread);
//            double y_new = minDisplay + (((xydata[i].getY() - min_y) / (max_y - min_y)) * 0.99 * displaySpread);
//            double z_new = minDisplay + (((xydata[i].getZ() - min_z) / (max_z - min_z)) * 0.99 * displaySpread);
//            float[] f = new float[3];
//            f[0] = (float) x_new;
//            f[1] = (float) y_new;
//            f[2] = (float) z_new;
//            xydata[i] = new Point3f(f);
//            setLocation(i, xydata[i]);
//        }
//    }

//    /* (non-Javadoc)
//     * @see edu.uci.ics.jung.visualization.layout.AbstractLayout#setSize(java.awt.Dimension)
//     */
//    @Override
//    public void setSize(Dimension size) {
//        setInitializer(new RandomLocationTransformer<Integer>(size));
//        super.setSize(size);
//    }
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
