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
public class WeightedMajorizationClustered {

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
    int[][] clIds;
    int[] numCls;
    Random r;
    static int maxIteration = 5;

    public WeightedMajorizationClustered() {
    }

    public void setGroundTruth(double[][] groundTruth) {
        this.groundTruth = groundTruth;
    }

    public WeightedMajorizationClustered(Graph g, int d, int randomSeed, int[][] clId, int[] numCl) {
        this.g = g;
        this.d = d;
        this.clIds = clId;
        this.numCls = numCl;
        this.randomSeed = randomSeed;
        this.maxSigma = maxSigma;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        currentIteration = 0;
        coord = new double[n][d];
        r = new Random(randomSeed);
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

    public WeightedMajorizationClustered(Graph g, int d, int randomSeed, double maxSigma) {
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

    public WeightedMajorizationClustered(Graph g, int d, double maxSigma, double[][] groundTruth) {
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
        int maxIter = 10000;
        for (int i = 0; i < maxIter; i++) {
            update();
            if (i % 200 == 0) {
                DataUtils du = new DataUtils();
                String fn = "c_" + i; 
                du.saveAsMatlab(coord, "coord", "result.mat");
            }
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

        int index = r.nextInt(clIds.length);

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
//        if (s.sigma > maxSigma) {
//            s.sigma = maxSigma;
//        }
        sigma = s.sigma;

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
                double[] dw = s.parabola(x, g.isNeighbor(i, j));
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
        if (currentIteration % 10 == 0) {
            System.out.println(currentIteration + " " + s.sigma + " " + aktCost + " " + index);
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

        lastCost = aktCost;

        //copy coord, distance and weights
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < d; j++) {
                coord[i][j] = coord_new[i][j];
            }
        }

//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < i; j++) {
//                weights[i][j] = weights[j][i] = weights_new[i][j];
//                distances[i][j] = distances[j][i] = dist_new[i][j];
//            }
//        }
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

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
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

}
