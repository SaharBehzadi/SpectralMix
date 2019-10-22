/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import java.util.ConcurrentModificationException;
import java.util.List;
import java.util.Random;
import mdsj.MDSJ;
import mdsj.StressMinimization;

/**
 *
 * @author claudia.plant
 */
public class WeightedMajorization {

    int currentIteration;
    int bestIteration;
    int randomSeed;
    double[][] coord;
    double[][] distances;
    double[][] weights;
    double[][] groundTruth; //ground truth coordinates
    int n; //numObj
    int d; //embedding dim
    int m;
    Graph g;
    boolean isoInit = false;
    boolean groundTruthInit = false;
    boolean randomInit = true;
    boolean mdsInit = false;
    GraphCompression gk;
    double costWithoutEmbedding;
    double bestCost;
    double paramCost;
    double savedBits;
    double lastCost;
    double[][] bestCoord;
    double sigma;
    double maxSigma;
    static int maxIteration = 5;

    public WeightedMajorization() {
    }

    public void setGroundTruth(double[][] groundTruth) {
        this.groundTruth = groundTruth;
    }

    public WeightedMajorization(Graph g, int d, int randomSeed, double maxSigma) {
        this.g = g;
        this.d = d;
        this.randomSeed = randomSeed;
        this.maxSigma = maxSigma;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        currentIteration = 0;
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
            //Random r = new Random();
            Random r = new Random(randomSeed);
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
        bestCoord = coord.clone();
        weights = new double[n][n];
        gk = new GraphCompression(g, coord);
        costWithoutEmbedding = gk.codingCostNoEmbedding();
        //System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        gk.setCoord(coord);
        bestCost = gk.mdlFunctionSimpleSigmoid();
        paramCost = gk.paramCostBic();
        //System.out.println("pCostsBic: " + paramCost);
        savedBits = costWithoutEmbedding - (bestCost + paramCost);
        //System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
    }

    public WeightedMajorization(Graph g, int d, double maxSigma, double[][] groundTruth) {
        this.g = g;
        this.d = d;
        this.randomSeed = randomSeed;
        this.maxSigma = maxSigma;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        currentIteration = 0;
        coord = new double[n][d];
        //init coord with isomap

        coord = groundTruth;
        distances = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
            }
        }

        bestCoord = coord.clone();
        weights = new double[n][n];
        gk = new GraphCompression(g, coord);
        costWithoutEmbedding = gk.codingCostNoEmbedding();
      //  System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        gk.setCoord(coord);
        bestCost = gk.mdlFunctionSimpleSigmoid();
        paramCost = gk.paramCostBic();
       // System.out.println("pCostsBic: " + paramCost);
        savedBits = costWithoutEmbedding - (bestCost + paramCost);
        //System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
    }

    public double[][] run() {
        int maxIter = 2000;
        for (int i = 0; i < maxIter; i++) {
            update();
        }
        return bestCoord;

    }

    //8.5.15: for hierarchical embedding: init with coord run until convergence
    public void run(double convConst) {

        boolean converged = false;
        double aktCost = Double.MAX_VALUE;
        while ((!converged || currentIteration < 5) && currentIteration < maxIteration) {
            //while (!converged) {
            //DEBUG
//            if(currentIteration == 173)
//                System.out.println("m");
            //DEBUG
            update();
            // System.out.println("converged after: " + currentIteration + " " + sigma + " " + lastCost);
            if ((aktCost - bestCost) > convConst) {
                aktCost = bestCost;
            } else {
                converged = true;
            }
        }
        //System.out.println(currentIteration + " iterations. Sigma: " + sigma + " cost: " + bestCost + " saved bits: " + savedBits + " entropy: " + costWithoutEmbedding + " parameterCosts: " + paramCost);
    }

    //9.1.15
    public void writeSortedAdjacenyMatrix() {
    }

    private void update() {

        double minWeight = 1E-9;
        double[][] weights_new = new double[n][n];
        double[][] dist_new = new double[n][n];

        //compute sigmoid
        PairOwn[] p = new PairOwn[m];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), g.isNeighbor(i, j));
            }
        }
//        double meanDist = 0.0;
//        for (int i = 0; i < p.length; i++) {
//            meanDist += p[i].dist;
//        }
//        meanDist /= (double) p.length;
//        System.out.println("meanDist " + meanDist);
        SimpleSigmoid s = new SimpleSigmoid(p);
//        if (s.sigma > 5.0) {
//            s.sigma = 5.0;
//        }
//            if (loop > ) {
////                if(loop == 22)
////                    System.out.println("m");
//                s = new SimpleSigmoid(1.0);
//            }
//            if (loop > 50) {
//                s = new SimpleSigmoid(p);
//            }

        double aktCost = s.costAllPairs(p);
        if (s.sigma > maxSigma) {
            s.sigma = maxSigma;
        }
        sigma = s.sigma;

        //System.out.println(currentIteration + " " + s.mu + " " + s.sigma + " Cost=" + aktCost);
        // if(cheatVar)
        //  s.sigma = Math.max(s.sigma, 5.0);
        // System.out.println("variance before update: " + s.sigma);
        //compute new weights and dists by setting derivatives of MDS error function and cost function equal
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

                }

            }
        }
        //Test
        //  weights_new = scale(weights_new);
        //dist_new = scale(dist_new);
        DataUtils du = new DataUtils();

        String name = "results_" + currentIteration;

        String name1 = "dist_" + currentIteration;
        String name2 = "weights_" + currentIteration;
        IO ea = new IO();
//        ea.writeDoubleToMatlab(dist_new, name1);
//        ea.writeDoubleToMatlab(weights_new, name2);
        //Test

        //do embedding
//        StressMinimization sm = new StressMinimization(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
//        String ss = sm.iterate(1);
//        //System.out.println(ss);
//        double[][] bla = sm.getPositions();
        WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dist_new, new Matrix(coord).transpose().getArrayCopy(), weights_new);
        wp.iterate();
        double[][] bla = wp.getPositions();

        //double[][] bla = MDSJ.stressMinimization(dist_new, weights_new, 2);
        double[][] coord_new = new Matrix(bla).transpose().getArrayCopy();
        gk.setCoord(coord_new);
        //double aktCost = gk.mdlFunctionSimpleSigmoid(false);

//        String name3 = "coord_" + currentIteration;
//        ea.writeDoubleToMatlab(coord, name3);
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

        if (aktCost < bestCost) {
            bestCost = aktCost;
            savedBits = costWithoutEmbedding - (aktCost + paramCost);
            bestIteration = currentIteration;
            //System.out.println(loop + ", " + aktCost + " saved bits: " + savedBits);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    bestCoord[i][j] = coord[i][j];
                }
                //xydata[i].setLocation(coord[i][0], coord[i]e[1]);
            }
        }

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

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }
}
