/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import java.util.Collection;
import java.util.ConcurrentModificationException;
import java.util.List;
import java.util.Random;
import mdsj.MDSJ;
import mdsj.StressMinimization;

/**
 *
 *
 */
public class Grid {

    private int currentIteration;
    private double bestCost; //overall best cost
    private double lastCost = Double.MAX_VALUE; //cost in previous interation
    private int bestIteration; //iteration with the best cost
    private int d = 2;
    private int[] vertices;
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
    boolean isoInit = true;
    boolean mdsInit = false;
    boolean groundTruthInit = false;
    boolean randomInit = false;
    int convergeCounter = 0;
    boolean varianceFree = true;
    Graph g;
    GraphCompression gk;
    int n;
    int m;
    int numBins;
    double maxSigma;
    boolean verbose = true;
    static int maxIteration = 10000;

    public Grid(Graph g) {
        this.g = g;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        //this.maxSigma = maxSigma;
    }

    public Grid(Graph g, double maxSigma) {
        this.g = g;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        this.maxSigma = maxSigma;
    }

    public void setTau(double tau) {
        this.tau = tau;
    }

    //8.5.15: for hierarchical embedding: init with coord run until convergence
    public void run(double[][] coord, double convConst) {
        this.coord = coord;
        init();

        boolean converged = false;
        double aktCost = Double.MAX_VALUE;
        int minIter = 10000;
        while ((!converged || currentIteration < minIter) && currentIteration < maxIteration) {
            //while (!converged) {
            //DEB
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
            if (currentIteration % 100 == 0) {
               System.out.println(currentIteration + " " + aktCost + " " + sigma);
                DataUtils du = new DataUtils();
                int rCost = (int) (Math.floor(aktCost));
                String fName = Integer.toString(currentIteration) + "_" + Integer.toString(rCost) +  ".mat";
                du.saveAsMatlab(coord, "coord", fName);
            }
//            if (currentIteration % 1000 == 0) {
//                DataUtils du = new DataUtils();
//                int rCost = (int) Math.floor(aktCost);
//                String fName = Integer.toString(currentIteration) + "_" + Integer.toString(rCost) + "_" + sigma + ".mat";
//                du.saveAsMatlab(coord, "coord", fName);
//            }
        }
        if (verbose) {

            System.out.println(currentIteration + " iterations. Sigma: " + sigma + " cost: " + bestCost + " saved bits: " + savedBits + " entropy: " + costWithoutEmbedding + " parameterCosts: " + paramCost);
        }
    }
    
    //8.5.15: for hierarchical embedding: init with coord run until convergence
    public void run(double convConst) {
        initialize();

        boolean converged = false;
        double aktCost = Double.MAX_VALUE;
        int minIter = 50;
        while ((!converged || currentIteration < minIter) && currentIteration < maxIteration) {
            //while (!converged) {
            //DEB
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
//            if (currentIteration % 100 == 0) {
//               System.out.println(currentIteration + " " + aktCost + " " + sigma);
//                DataUtils du = new DataUtils();
//                int rCost = (int) (Math.floor(aktCost));
//                String fName = Integer.toString(currentIteration) + "_" + Integer.toString(rCost) +  ".mat";
//                du.saveAsMatlab(coord, "coord", fName);
//            }
            if (currentIteration % 1000 == 0) {
//                DataUtils du = new DataUtils();
//                int rCost = (int) Math.floor(aktCost);
//                String fName = Integer.toString(currentIteration) + "_" + Integer.toString(rCost) + "_" + sigma + ".mat";
//                du.saveAsMatlab(coord, "coord", fName);
            }
        }
        if (verbose) {

            System.out.println(currentIteration + " iterations. Sigma: " + sigma + " cost: " + bestCost + " saved bits: " + savedBits + " entropy: " + costWithoutEmbedding + " parameterCosts: " + paramCost);
        }
    }

    public void run() {
        DataUtils du = new DataUtils();
        double[][] costSigma = new double[maxIteration + 1][2];
        initialize();
        costSigma[currentIteration][0] = bestCost;
        costSigma[currentIteration][1] = sigma;
        System.out.println(currentIteration + " " + sigma + " " + bestCost);
        while (currentIteration < maxIteration) {
            update();
            costSigma[currentIteration][0] = lastCost;
            costSigma[currentIteration][1] = sigma;
            if (currentIteration % 100 == 0) {
                String name = "coord_" + currentIteration;
                du.saveAsMatlab(coord, "coord", name + ".mat");
                System.out.println(currentIteration + " " + sigma + " " + lastCost);
                // du.saveAsMatlab3(dist_new, weights_new, coord_new, "dist", "weight", "coord", name+".mat");
                // ea.writeDoubleToMatlab(coord, name3);

            }

        }
        du.saveAsMatlab(costSigma, "costSigma", "costSigma.mat");
    }

    public double[][] getCoord() {
        return coord;
    }

    public void init() {
        currentIteration = 0;
        vertices = new int[n];
        for (int i = 0; i < n; i++) {
            vertices[i] = i;
        }
        positions = new double[n][d];
        gk = new GraphCompression(g, coord);
        costWithoutEmbedding = gk.codingCostNoEmbedding();
        //System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);

        bestDb = new double[n][d];
        for (int i = 0; i < n; i++) {
            bestDb[i][0] = coord[i][0];
            bestDb[i][1] = coord[i][1];
        }
        gk.setCoord(coord);
        // = 0;
        bestCost = gk.mdlFunctionSimpleSigmoid();
        paramCost = gk.paramCostBic();
        mu = gk.s.mu;
        sigma = gk.s.sigma;
        //tau = 1E-3;  //1E-6
        tau = 1E-4;
        bestIteration = 0;
        //System.out.println("pCostsBic: " + paramCost);
        savedBits = costWithoutEmbedding - (bestCost + paramCost);
        //System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
    }

    public void initialize() {
        currentIteration = 0;
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

        bestDb = new double[n][d];
        for (int i = 0; i < n; i++) {
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
        bestIteration = 0;
        System.out.println("pCostsBic: " + paramCost);
        savedBits = costWithoutEmbedding - (bestCost + paramCost);
        System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
    }

    public double getBestCost() {
        return bestCost;
    }

    public double[][] getBestDb() {
        return bestDb;
    }

    public int getBestIteration() {
        return bestIteration;
    }

    public double getSavedBits() {
        return savedBits;
    }

    private void update() {
        if(currentIteration == 25)
            System.out.println("m");
        DataUtils du = new DataUtils();

        double cutOff = gridThreshold();
        //cutOff = 100;
        double[][] minMax = du.minMax(coord);
        int numBinsX = (int) Math.ceil((minMax[0][1] - minMax[0][0]) / cutOff);
        int numBinsY = (int) Math.ceil((minMax[1][1] - minMax[1][0]) / cutOff);
        // System.out.println(cutOff + " " + numBinsX + " " + numBinsY);
//        numBinsX = Math.min(numBinsX, 1);
//        numBinsY = Math.min(numBinsY, 1);

        int[][] numPoints = new int[numBinsX][numBinsY];
        for (int i = 0; i < n; i++) {
            //compute grid cell of point
            //first bin: [0, 0.01[, last bin [0.9, 1.0]
//            if(i == 186)
//                System.out.println("m");
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

        //System.out.println("numBinsX " + numBinsX + " numBinsY " + numBinsY);
        //compute sigmoid
        PairOwn[] p = new PairOwn[m];
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

        SimpleSigmoid s = new SimpleSigmoid(p);
        if (s.sigma > maxSigma) {
            s.sigma = maxSigma;
        }
        mu = s.mu;
        sigma = s.sigma;
        //Frisieren

        double aktCost = s.costAllPairs(p);
        //System.out.println(aktCost);

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
        System.out.println(currentIteration + ", " + aktCost + " " + s.sigma);
        lastCost = aktCost;

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

        if (aktCost < bestCost) {
            bestCost = aktCost;
            bestIteration = currentIteration;
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
    }

    private double gridThreshold() {
        if (Math.pow(tau, 2) * Math.pow(sigma, 4) > 0.05) {
            return 1.33;
            
        }
        double xiOld = 1.0;
        double xiNew = 0.0;
        boolean converged = false;
        double convConst = 1E-3;
        double xi = 1.0;
        int iter = 0;
        while (!converged) {
            iter++;
            double xiTerm = Math.pow(tau, 2) / Math.pow(xiOld, 2);
            double product = 4 * Math.PI * xiTerm * Math.pow(Math.log(2), 2) * Math.pow(sigma, 4);
            xiNew = Math.sqrt(-Math.log(product));
            if (Math.abs(xiOld - xiNew) < convConst || iter > 2000) {
                converged = true;
            }
            xiOld = xiNew;
        }

        double dist = Math.sqrt(2.0) * sigma * xiNew + mu;
        return dist;
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
