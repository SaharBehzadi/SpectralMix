/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.locks.ReentrantLock;

/**
 *
 * @author plantc59cs
 */
public class IncrementalParallelGrid {

    private double bestCost; //overall best cost
    private double lastCost = Double.MAX_VALUE; //cost in previous interation
    public int d = 2;
    public VertexLocks[] vertices;
    double costWithoutEmbedding;
    double[][] bestCoord;
//    double[][] distances;
    double[][] coord;
    // double[][] weights;
    double[][] groundTruth; //ground truth coordinates
    boolean randomInit = true;
    int convergeCounter = 0;
    boolean varianceFree = false;
    boolean useGrid;
    Graph g;
    GraphCompression gk;
    SimpleSigmoid s;
    Random r;
    // int d;
    int n;
    int m;
    Dimension dd; //size of display
    //grid
    int[][][] pId;
    int numBinsX; //grid
    int numBinsY; //grid
    double[][] minMax; //grid
    double cutOff; //grid
    static double tau = 0.1;
    static double mu = 1.5;
    //new stuff
    boolean[] isEdge;
    int[] nonEdgeFirst;
    int[] nonEdgeSecond;
    double[][] nomEdge;
    double[][] nomEdgeR;
    double[] denomEdge;
    double[][] nomVertex;
    double[] denomVertex;
    ReentrantLock[] locks;
    boolean[] valid;
    boolean verbose = true;
    boolean writeResult = true;
    boolean display = false;

    public IncrementalParallelGrid(Graph<Integer, Integer> g) {
        this.g = g;
        r = new Random(1);
        n = g.getVertexCount();
        coord = new double[n][d];
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < d; j++) {
                coord[i][j] = r.nextDouble();
            }
        }
        setUpGrid();
        //collect which non-edges need to be updated
        Vector<Integer> first = new Vector<Integer>();
        Vector<Integer> second = new Vector<Integer>();
        for (int i = 0; i < numBinsX; i++) {
            for (int j = 0; j < numBinsY; j++) {
                for (int a = 0; a < pId[i][j].length; a++) {
                    for (int b = 0; b < a; b++) {
                        if (!g.isNeighbor(pId[i][j][a], pId[i][j][b])) {
                            first.add(pId[i][j][a]);
                            second.add(pId[i][j][b]);
                        }
                    }
                }
                //pId[i + 1][j][b]
                if (i + 1 < numBinsX) {
                    for (int a = 0; a < pId[i][j].length; a++) {
                        for (int b = 0; b < pId[i + 1][j].length; b++) {
                            if (!g.isNeighbor(pId[i][j][a], pId[i + 1][j][b])) {
                                first.add(pId[i][j][a]);
                                second.add(pId[i + 1][j][b]);
                            }
                        }
                    }
                }
                //pId[i2][j + 1][b]
                if (j + 1 < numBinsY) {
                    for (int i2 = Math.max(0, i - 1); i2 < Math.min(numBinsX, i + 2); i2++) {
                        for (int a = 0; a < pId[i][j].length; a++) {
                            for (int b = 0; b < pId[i2][j + 1].length; b++) {
                                if (!g.isNeighbor(pId[i][j][a], pId[i2][j + 1][b])) {
                                    first.add(pId[i][j][a]);
                                    second.add(pId[i2][j + 1][b]);
                                }
                            }
                        }
                    }
                }
            }
        }

        m = g.getEdgeCount() + nonEdgeFirst.length;
        nonEdgeFirst = new int[first.size()];
        nonEdgeSecond = new int[second.size()];
        for (int i = 0; i < nonEdgeFirst.length; i++) {
            nonEdgeFirst[i] = first.elementAt(i);
            nonEdgeSecond[i] = second.elementAt(i);
        }
        isEdge = new boolean[m];
        for (int i = 0; i < g.getEdgeCount(); i++) {
            isEdge[i] = true;
        }
        valid = new boolean[m];
        nomEdge = new double[m][d];
        nomEdgeR = new double[m][d];
        denomEdge = new double[m];
        nomVertex = new double[n][d];
        denomVertex = new double[n];
        locks = new ReentrantLock[n];
        for (int i = 0; i < n; i++) {
            locks[i] = new ReentrantLock(true);
        }

        gk = new GraphCompression(g, coord);
        s = new SimpleSigmoid(1.0);
        gk.setCoord(coord);
//        loop = 0;
        bestCoord = new double[n][d];
        bestCost = gk.mdlFunctionSimpleSigmoid();
        if (display) {
            Visualization v = new Visualization(g);
            v.displayCoordNew(coord, "init");
        }
        if (verbose) {
            System.out.println("cost afte init: " + bestCost + " " + gk.s.sigma);
        }
    }

    public IncrementalParallelGrid(Graph<Integer, Integer> g, double[][] cc) {
        this.g = g;
        this.coord = cc;
        r = new Random(1);
        n = g.getVertexCount();
        setUpGrid();
        //collect which non-edges need to be updated
        Vector<Integer> first = new Vector<Integer>();
        Vector<Integer> second = new Vector<Integer>();
        for (int i = 0; i < numBinsX; i++) {
            for (int j = 0; j < numBinsY; j++) {
                for (int a = 0; a < pId[i][j].length; a++) {
                    for (int b = 0; b < a; b++) {
                        if (!g.isNeighbor(pId[i][j][a], pId[i][j][b])) {
                            first.add(pId[i][j][a]);
                            second.add(pId[i][j][b]);
                        }
                    }
                }
                //pId[i + 1][j][b]
                if (i + 1 < numBinsX) {
                    for (int a = 0; a < pId[i][j].length; a++) {
                        for (int b = 0; b < pId[i + 1][j].length; b++) {
                            if (!g.isNeighbor(pId[i][j][a], pId[i + 1][j][b])) {
                                first.add(pId[i][j][a]);
                                second.add(pId[i + 1][j][b]);
                            }
                        }
                    }
                }
                //pId[i2][j + 1][b]
                if (j + 1 < numBinsY) {
                    for (int i2 = Math.max(0, i - 1); i2 < Math.min(numBinsX, i + 2); i2++) {
                        for (int a = 0; a < pId[i][j].length; a++) {
                            for (int b = 0; b < pId[i2][j + 1].length; b++) {
                                if (!g.isNeighbor(pId[i][j][a], pId[i2][j + 1][b])) {
                                    first.add(pId[i][j][a]);
                                    second.add(pId[i2][j + 1][b]);
                                }
                            }
                        }
                    }
                }
            }
        }

        nonEdgeFirst = new int[first.size()];
        nonEdgeSecond = new int[second.size()];
        for (int i = 0; i < nonEdgeFirst.length; i++) {
            nonEdgeFirst[i] = first.elementAt(i);
            nonEdgeSecond[i] = second.elementAt(i);
        }
        m = g.getEdgeCount() + nonEdgeFirst.length;
        isEdge = new boolean[m];
        for (int i = 0; i < g.getEdgeCount(); i++) {
            isEdge[i] = true;
        }
        valid = new boolean[m];
        nomEdge = new double[m][d];
        nomEdgeR = new double[m][d];
        denomEdge = new double[m];
        nomVertex = new double[n][d];
        denomVertex = new double[n];
        locks = new ReentrantLock[n];
        for (int i = 0; i < n; i++) {
            locks[i] = new ReentrantLock(true);
        }
        gk = new GraphCompression(g, coord);

        gk.setCoord(coord);
//        loop = 0;
        bestCoord = new double[n][d];
        bestCost = gk.mdlFunctionSimpleSigmoid();
        if (display) {
            Visualization v = new Visualization(g);
            v.displayCoordNew(coord, "init");
        }
        if (verbose) {
            System.out.println("cost afte init: " + bestCost + " " + gk.s.sigma);
        }
    }

    private void computeSigmoid() {
        int numElt = n * (n - 1) / 2;
        PairOwn[] p = new PairOwn[numElt];
        //System.out.println(numElt);
        for (int i = 1; i < n; i++) { //why i = 1?
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), g.isNeighbor(i, j));
            }
        }

        s = new SimpleSigmoid(p);
    }

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }

    private void collectRelevantNotEdges() {
        //collect which non-edges need to be updated
        Vector<Integer> first = new Vector<Integer>();
        Vector<Integer> second = new Vector<Integer>();
        for (int i = 0; i < numBinsX; i++) {
            for (int j = 0; j < numBinsY; j++) {
                for (int a = 0; a < pId[i][j].length; a++) {
                    for (int b = 0; b < a; b++) {
                        if (!g.isNeighbor(pId[i][j][a], pId[i][j][b])) {
                            first.add(pId[i][j][a]);
                            second.add(pId[i][j][b]);
                        }
                    }
                }
                //pId[i + 1][j][b]
                if (i + 1 < numBinsX) {
                    for (int a = 0; a < pId[i][j].length; a++) {
                        for (int b = 0; b < pId[i + 1][j].length; b++) {
                            if (!g.isNeighbor(pId[i][j][a], pId[i + 1][j][b])) {
                                first.add(pId[i][j][a]);
                                second.add(pId[i + 1][j][b]);
                            }
                        }
                    }
                }
                //pId[i2][j + 1][b]
                if (j + 1 < numBinsY) {
                    for (int i2 = Math.max(0, i - 1); i2 < Math.min(numBinsX, i + 2); i2++) {
                        for (int a = 0; a < pId[i][j].length; a++) {
                            for (int b = 0; b < pId[i2][j + 1].length; b++) {
                                if (!g.isNeighbor(pId[i][j][a], pId[i2][j + 1][b])) {
                                    first.add(pId[i][j][a]);
                                    second.add(pId[i2][j + 1][b]);
                                }
                            }
                        }
                    }
                }
            }
        }

        nonEdgeFirst = new int[first.size()];
        nonEdgeSecond = new int[second.size()];
        for (int i = 0; i < nonEdgeFirst.length; i++) {
            nonEdgeFirst[i] = first.elementAt(i);
            nonEdgeSecond[i] = second.elementAt(i);
        }
        m = g.getEdgeCount() + nonEdgeFirst.length;
        isEdge = new boolean[m];
        for (int i = 0; i < g.getEdgeCount(); i++) {
            isEdge[i] = true;
        }
        valid = new boolean[m];
        nomEdge = new double[m][d];
        nomEdgeR = new double[m][d];
        denomEdge = new double[m];
    }

    private void setUpGrid() {
        computeSigmoid();
        DataUtils du = new DataUtils();
        cutOff = gridThreshold();
        //cutOff = 100;
        minMax = du.minMax(coord);
        numBinsX = (int) Math.ceil((minMax[0][1] - minMax[0][0]) / cutOff);
        numBinsY = (int) Math.ceil((minMax[1][1] - minMax[1][0]) / cutOff);
        //debug
//        numBinsX = 1;
//        numBinsY = 1;
        useGrid = numBinsX > 1 && numBinsY > 1;

        //System.out.println(cutOff + " " + numBinsX + " " + numBinsY);
//        numBinsX = Math.min(numBinsX, 1);
//        numBinsY = Math.min(numBinsY, 1);
        int[][] numPoints = new int[numBinsX][numBinsY];
        for (int i = 0; i < n; i++) {
            //compute grid cell of point
            //first bin: [0, 0.01[, last bin [0.9, 1.0]
            int iIndex = Math.min((int) Math.floor((coord[i][0] - minMax[0][0]) / cutOff), numBinsX - 1);
            int jIndex = Math.min((int) Math.floor((coord[i][1] - minMax[1][0]) / cutOff), numBinsY - 1);
//              
            numPoints[iIndex][jIndex]++;
        }

        pId = new int[numBinsX][numBinsY][];
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
    }

    private double gridThreshold() {
        if (Math.pow(tau, 2) * Math.pow(s.sigma, 4) > 0.05) {
            return 1.33;

        }
        double xiOld = 1.0;
        double xiNew = 0.0;
        boolean converged = false;
        double cc = 1E-3;
        double xi = 1.0;
        int ii = 0;
        while (!converged) {
            ii++;
            double xiTerm = Math.pow(tau, 2) / Math.pow(xiOld, 2);
            double product = 4 * Math.PI * xiTerm * Math.pow(Math.log(2), 2) * Math.pow(s.sigma, 4);
            xiNew = Math.sqrt(-Math.log(product));
            if (Math.abs(xiOld - xiNew) < cc || ii > 2000) {
                converged = true;
            }
            xiOld = xiNew;
        }

        double dist = Math.sqrt(2.0) * s.sigma * xiNew + mu;
        return dist;
    }

    public void lock4o(int la, int lb, int lc, int ld) {
        // System.out.println(la+" "+lb+" "+lc+" "+ld);
        locks[la].lock(); // BITTE ANPASSEN !!!
        locks[lb].lock();
        locks[lc].lock();
        locks[ld].lock();
        //System.out.println(".");
    }

    public void lock4u(int la, int lb, int lc, int ld) {
        // System.out.println(la+" "+lb+" "+lc+" "+ld);
        if (la < lb) {
            if (lc < ld) {
                if (lb < lc) {
                    lock4o(la, lb, lc, ld);
                } else if (la < lc) {
                    if (lb < ld) {
                        lock4o(la, lc, lb, ld);
                    } else {
                        lock4o(la, lc, ld, lb);
                    }
                } else if (lb < ld) {
                    lock4o(lc, la, lb, ld);
                } else if (la < ld) {
                    lock4o(lc, la, ld, lb);
                } else {
                    lock4o(lc, ld, la, lb);
                }
            } else // vertausche Rollen von lc und ld
             if (lb < ld) {
                    lock4o(la, lb, ld, lc);
                } else if (la < ld) {
                    if (lb < lc) {
                        lock4o(la, ld, lb, lc);
                    } else {
                        lock4o(la, ld, lc, lb);
                    }
                } else if (lb < lc) {
                    lock4o(ld, la, lb, lc);
                } else if (la < lc) {
                    lock4o(ld, la, lc, lb);
                } else {
                    lock4o(ld, lc, la, lb);
                }
        } else // vertausche Rollen von la und lb
         if (lc < ld) {
                if (la < lc) {
                    lock4o(lb, la, lc, ld);
                } else if (lb < lc) {
                    if (la < ld) {
                        lock4o(lb, lc, la, ld);
                    } else {
                        lock4o(lb, lc, ld, la);
                    }
                } else if (la < ld) {
                    lock4o(lc, lb, la, ld);
                } else if (lb < ld) {
                    lock4o(lc, lb, ld, la);
                } else {
                    lock4o(lc, ld, lb, la);
                }
            } else // vertausche Rollen von lc und ld
             if (la < ld) {
                    lock4o(lb, la, ld, lc);
                } else if (lb < ld) {
                    if (la < lc) {
                        lock4o(lb, ld, la, lc);
                    } else {
                        lock4o(lb, ld, lc, la);
                    }
                } else if (la < lc) {
                    lock4o(ld, lb, la, lc);
                } else if (lb < lc) {
                    lock4o(ld, lb, lc, la);
                } else {
                    lock4o(ld, lc, lb, la);
                }
    }

    public void firstRound(int numT) {
        int numEdgesPerThread = m / numT;
        int remain = m - (numT * numEdgesPerThread);
        int counter = 0;
        int counterNonEdges = 0;
        ExecutorService executor = Executors.newFixedThreadPool(numT);
        Collection collection = new ArrayList();

        for (int i = 0; i < numT - 1; i++) {
            Pair[] edges = new Pair[numEdgesPerThread];
            int[] eIds = new int[numEdgesPerThread];
            for (int j = 0; j < edges.length; j++) {
                if (counter < g.getEdgeCount()) {
                    edges[j] = g.getEndpoints(counter);
                } else {
                    edges[j] = new Pair(nonEdgeFirst[counterNonEdges], nonEdgeSecond[counterNonEdges]);
                    counterNonEdges++;
                }
                eIds[j] = counter;
                counter++;
            }

            StartEdgeUpdateThreadGrid t = new StartEdgeUpdateThreadGrid(edges, eIds, this, r.nextInt());

            collection.add(t);
        }

        Pair[] edges = new Pair[numEdgesPerThread + remain];
        int[] eIds = new int[numEdgesPerThread + remain];
        for (int j = 0; j < edges.length; j++) {
            if (counter < g.getEdgeCount()) {
                edges[j] = g.getEndpoints(counter);
            } else {
                edges[j] = new Pair(nonEdgeFirst[counterNonEdges], nonEdgeSecond[counterNonEdges]);
                counterNonEdges++;
            }
            eIds[j] = counter;
            counter++;
        }
        //EdgeUpdateThread t = new EdgeUpdateThread(edges, eIds, this, iterPerT, r.nextInt());
        StartEdgeUpdateThreadGrid t = new StartEdgeUpdateThreadGrid(edges, eIds, this, r.nextInt());

        collection.add(t);

        try {
            executor.invokeAll(collection);
        } catch (Exception e) {
            e.printStackTrace();
        }

        executor.shutdown();

    }

    public void runUpdates(int numT, int numIterPerT) {
        int numEdgesPerThread = m / numT;
        int remain = m - (numT * numEdgesPerThread);
        int counter = 0;
        int counterNonEdges = 0;
        ExecutorService executor = Executors.newFixedThreadPool(numT);
        Collection collection = new ArrayList();

        for (int i = 0; i < numT - 1; i++) {
            Pair[] edges = new Pair[numEdgesPerThread];
            int[] eIds = new int[numEdgesPerThread];
            for (int j = 0; j < edges.length; j++) {
                if (counter < g.getEdgeCount()) {
                    edges[j] = g.getEndpoints(counter);
                } else {
                    edges[j] = new Pair(nonEdgeFirst[counterNonEdges], nonEdgeSecond[counterNonEdges]);
                    counterNonEdges++;
                }
                eIds[j] = counter;
                counter++;
            }
            EdgeUpdateThreadGrid t = new EdgeUpdateThreadGrid(edges, eIds, this, numIterPerT);
            collection.add(t);
        }
        Pair[] edges = new Pair[numEdgesPerThread + remain];
        int[] eIds = new int[numEdgesPerThread + remain];
        for (int j = 0; j < edges.length; j++) {
            if (counter < g.getEdgeCount()) {
                edges[j] = g.getEndpoints(counter);
            } else {
                edges[j] = new Pair(nonEdgeFirst[counterNonEdges], nonEdgeSecond[counterNonEdges]);
                counterNonEdges++;
            }
            eIds[j] = counter;
            counter++;
        }
        //EdgeUpdateThread t = new EdgeUpdateThread(edges, eIds, this, iterPerT, r.nextInt());
        EdgeUpdateThreadGrid t = new EdgeUpdateThreadGrid(edges, eIds, this, numIterPerT);

        collection.add(t);

        try {
            executor.invokeAll(collection);
        } catch (Exception e) {
            e.printStackTrace();
        }

        executor.shutdown();

    }

    private boolean checkLocks() {
        boolean allFree = true;
        for (int i = 0; i < locks.length; i++) {
            if (locks[i].isLocked() == true) {
                allFree = false;
            }
        }
        return allFree;
    }

    public void run(int iter, int numUpdatesPerEdge, int numT) {
        double startTime = System.currentTimeMillis();
        firstRound(numT);
        //boolean allFree = checkLocks();
        Visualization v = new Visualization(g);

        //int numUpdatesPerT = (n / numT) * numUpdatesPerEdge;
        //should be
        int numUpdatesPerT = (m / numT) * numUpdatesPerEdge;
        if (verbose) {
            System.out.println(m + " entries, corresponds to " + numUpdatesPerT + " updates per thread.");
        }
        for (int i = 0; i < iter; i++) {
            if (i > 0) {
                setUpGrid();
                collectRelevantNotEdges();
                firstRound(numT);
                numUpdatesPerT = (m / numT) * numUpdatesPerEdge;
                //System.out.println(numUpdatesPerT);
            }
            runUpdates(numT, numUpdatesPerT);
            if (verbose) {
                //if (i % 100 == 0) {          
                gk.setCoord(coord);
                System.out.println(i + " " + gk.mdlFunctionSimpleSigmoid() + " " + gk.s.sigma);
                //}
            }
           // if (writeResult && i %100 == 0) {
           if(writeResult){
                String s = new Integer(i).toString();
                DataUtils du = new DataUtils();
                du.saveAsMatlab(coord, "coord", "result" + s + ".mat");
           }
            //}
            if (display) {
                String s = new Integer(i).toString();
                //du.saveAsMatlab3(coord, nomVertex, denomV, "coord", "nomVertex", "denomVertex", s + ".mat");
                v.displayCoordNew(coord, s);
            }

        }
        double endTime = System.currentTimeMillis();
        if (verbose) {
            System.out.println("cost after refinement: " + gk.mdlFunctionSimpleSigmoid() + " " + gk.s.sigma);
            NumberFormat formatter = new DecimalFormat("#0.00000");
            System.out.println("runtime " + formatter.format((endTime - startTime) / 1000d) + " seconds");

        }

        if (display) {
            v.displayCoordNew(coord, "finished");

        }
        if (verbose) {
            gk.setCoord(coord);
            bestCost = gk.mdlFunctionSimpleSigmoid();
            System.out.println("cost: " + bestCost + " " + gk.s.sigma);
        }
        if (writeResult) {
            DataUtils du = new DataUtils();
            du.saveAsMatlab(coord, "coord", "result.mat");
        }

    }

//    public void runUpdates(int numT, int iterPerT) {
//        //public EdgeUpdateThread(Pair<Integer>[] edges, int[] eIds, IncrementalParallel p, int maxIter, int seed) {
//        int numEdgesPerThread = m / numT;
//        int remain = m - (numT * numEdgesPerThread);
//        int counter = 0;
//        ExecutorService executor = Executors.newFixedThreadPool(numT);
//        Collection collection = new ArrayList();
//
//        for (int i = 0; i < numT - 1; i++) {
//            Pair[] edges = new Pair[numEdgesPerThread];
//            int[] eIds = new int[numEdgesPerThread];
//            for (int j = 0; j < edges.length; j++) {
//                edges[j] = g.getEndpoints(counter);
//                eIds[j] = counter;
//                counter++;
//            }
//            EdgeUpdateThread t = new EdgeUpdateThread(edges, eIds, this, iterPerT, r.nextInt());
//            //StartEdgeUpdateThread t = new StartEdgeUpdateThread(edges, eIds, this, r.nextInt());
//            collection.add(t);
//        }
//
//        Pair[] edges = new Pair[numEdgesPerThread + remain];
//        int[] eIds = new int[numEdgesPerThread + remain];
//        for (int j = 0; j < edges.length; j++) {
//            edges[j] = g.getEndpoints(counter);
//            eIds[j] = counter;
//            counter++;
//        }
//        EdgeUpdateThread t = new EdgeUpdateThread(edges, eIds, this, iterPerT, r.nextInt());
//        //StartEdgeUpdateThread t = new StartEdgeUpdateThread(edges, eIds, this, r.nextInt());
//
//        collection.add(t);
//
//        try {
//            executor.invokeAll(collection);
//            //System.out.println("finished");
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//
//        executor.shutdown();
//
//    }
}
