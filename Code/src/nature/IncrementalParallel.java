/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import attributedEmbedding.VisuAttributes;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.locks.ReentrantLock;

/**
 *
 * @author plantc59cs
 */
public class IncrementalParallel {

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
    Graph g;
    GraphCompression gk;
    SimpleSigmoid s;
    Random r;
    // int d;
    int n;
    int m;
    Dimension dd; //size of display
    //new stuff
    double[][] nomEdge;
    double[][] nomEdgeR;
    double[] denomEdge;
    int[] vNonEdge1;
    int[] vNonEdge2;
    double[][] nomNonEdge1;
    double[][] nomNonEdge1R;
    double[] denomNonEdge1;
    double[][] nomNonEdge2;
    double[][] nomNonEdge2R;
    double[] denomNonEdge2;
    double[][] nomVertex;
    double[] denomVertex;
    //required for refinement with clustering/product quantization
    //represent the notlinks which are updated in the next step
    int[][] notLinks; //store for each point numNotLink notLinks which are close by
    int numNotLinks;
    double[] avgDist; //the average distance to those points
    int[] activeNotLinks; //the current number of entries of the notLinks array
    //represent the clusters
    int k; //number of clusters; can be used to determine the locality
    //double[][] sum; //sum of the coordinates of the points; k x d
    //HashSet<Integer>[] clusterMembers; //for fast sampling of points in the same cluster
    Cluster[] clusters;
    int[] clusterID;
    //end of stuff required for refinement
    ReentrantLock[] locks;
    boolean[] valid;
    boolean verbose = true;
    boolean writeResult = true;
    boolean display = false;
    boolean computeCost = true;

    public IncrementalParallel(Graph<Integer, Integer> g) {
        this.g = g;
        this.k = k;
        r = new Random(1);
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        valid = new boolean[g.getEdgeCount()];
        nomEdge = new double[g.getEdgeCount()][d];
        nomEdgeR = new double[g.getEdgeCount()][d];
        denomEdge = new double[g.getEdgeCount()];
        vNonEdge1 = new int[g.getEdgeCount()];
        vNonEdge2 = new int[g.getEdgeCount()];
        nomNonEdge1 = new double[g.getEdgeCount()][d];
        nomNonEdge1R = new double[g.getEdgeCount()][d];
        denomNonEdge1 = new double[g.getEdgeCount()];
        nomNonEdge2 = new double[g.getEdgeCount()][d];
        nomNonEdge2R = new double[g.getEdgeCount()][d];
        denomNonEdge2 = new double[g.getEdgeCount()];
        nomVertex = new double[n][d];
        denomVertex = new double[n];
        //currentIteration = 0;
//        bestloop = -1;
//        loop = 0;
//        vertices = new VertexLocks[n];
//        for (int i = 0; i < n; i++) {
//            vertices[i] = new VertexLocks(i);
//        }
        locks = new ReentrantLock[n];
        for (int i = 0; i < n; i++) {
            locks[i] = new ReentrantLock(true);
        }
        coord = new double[n][d];
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < d; j++) {
                coord[i][j] = r.nextDouble();
            }
        }
        s = new SimpleSigmoid(1.0);
        bestCoord = new double[n][d];

        if (computeCost) {
            gk = new GraphCompression(g, coord);
            gk.setCoord(coord);
            bestCost = gk.mdlFunctionSimpleSigmoid();
        }
        if (display) {
            Visualization v = new Visualization(g);
            v.displayCoordNew(coord, "init");
        }
        if (verbose && computeCost) {
            System.out.println("cost afte init: " + bestCost + " " + gk.s.sigma);

        }
    }

    public void initNotLinks(int size) {
        this.numNotLinks = size;
        avgDist = new double[n];
        activeNotLinks = new int[n];
        notLinks = new int[n][numNotLinks];
        for (int i = 0; i < n; i++) {
            notLinks[i] = new int[numNotLinks];
        }
    }

    private void reInit(int seed) {
        r = new Random(seed);
        valid = new boolean[g.getEdgeCount()];
        nomEdge = new double[g.getEdgeCount()][d];
        nomEdgeR = new double[g.getEdgeCount()][d];
        denomEdge = new double[g.getEdgeCount()];
        vNonEdge1 = new int[g.getEdgeCount()];
        vNonEdge2 = new int[g.getEdgeCount()];
        nomNonEdge1 = new double[g.getEdgeCount()][d];
        nomNonEdge1R = new double[g.getEdgeCount()][d];
        denomNonEdge1 = new double[g.getEdgeCount()];
        nomNonEdge2 = new double[g.getEdgeCount()][d];
        nomNonEdge2R = new double[g.getEdgeCount()][d];
        denomNonEdge2 = new double[g.getEdgeCount()];
        nomVertex = new double[n][d];
        denomVertex = new double[n];
        //currentIteration = 0;
//        bestloop = -1;
//        loop = 0;
//        vertices = new VertexLocks[n];
//        for (int i = 0; i < n; i++) {
//            vertices[i] = new VertexLocks(i);
//        }
        locks = new ReentrantLock[n];
        for (int i = 0; i < n; i++) {
            locks[i] = new ReentrantLock(true);
        }
        coord = new double[n][d];
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < d; j++) {
                coord[i][j] = r.nextDouble();
            }
        }
        gk = new GraphCompression(g, coord);
        s = new SimpleSigmoid(1.0);
        gk.setCoord(coord);
//        loop = 0;
        bestCost = gk.mdlFunctionSimpleSigmoid();

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
            {
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
            }
        } else // vertausche Rollen von la und lb
        {
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
            {
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
        }
    }

    public void firstRound(int numT) {
        int numEdgesPerThread = g.getEdgeCount() / numT;
        int remain = g.getEdgeCount() - (numT * numEdgesPerThread);
        int counter = 0;
        ExecutorService executor = Executors.newFixedThreadPool(numT);
        Collection collection = new ArrayList();

        for (int i = 0; i < numT - 1; i++) {
            Pair[] edges = new Pair[numEdgesPerThread];
            int[] eIds = new int[numEdgesPerThread];
            for (int j = 0; j < edges.length; j++) {
                edges[j] = g.getEndpoints(counter);
                eIds[j] = counter;
                counter++;
            }
            //EdgeUpdateThread t = new EdgeUpdateThread(edges, eIds, this, iterPerT, r.nextInt());
            StartEdgeUpdateThread t = new StartEdgeUpdateThread(edges, eIds, this, r.nextInt());
            collection.add(t);
        }

        Pair[] edges = new Pair[numEdgesPerThread + remain];
        int[] eIds = new int[numEdgesPerThread + remain];
        for (int j = 0; j < edges.length; j++) {
            edges[j] = g.getEndpoints(counter);
            eIds[j] = counter;
            counter++;
        }
        //EdgeUpdateThread t = new EdgeUpdateThread(edges, eIds, this, iterPerT, r.nextInt());
        StartEdgeUpdateThread t = new StartEdgeUpdateThread(edges, eIds, this, r.nextInt());

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

    public double[][] getBestCoord() {
        return bestCoord;
    }

    public void refine(int numUpdatesPerEdge, int numT, int numClusters) {
        //initNotLinks(numNotLinks);
        int numTry = 5;
        int iter = 10;
        initClusters(numClusters, numTry, iter);
        double startTime = System.currentTimeMillis();
        int numUpdatesPerT = (g.getEdgeCount() / numT) * numUpdatesPerEdge;
        if (verbose) {
            System.out.println(numUpdatesPerT + " updates per thread.");
        }
        int numEdgesPerThread = g.getEdgeCount() / numT;
        int remain = g.getEdgeCount() - (numT * numEdgesPerThread);
        int counter = 0;
        ExecutorService executor = Executors.newFixedThreadPool(numT);
        Collection collection = new ArrayList();

        for (int i = 0; i < numT - 1; i++) {
            Pair[] edges = new Pair[numEdgesPerThread];
            int[] eIds = new int[numEdgesPerThread];
            for (int j = 0; j < edges.length; j++) {
                edges[j] = g.getEndpoints(counter);
                eIds[j] = counter;
                counter++;
            }
            EdgeRefinementThread t = new EdgeRefinementThread(edges, eIds, this, numUpdatesPerT, r.nextInt());
            //StartEdgeUpdateThread t = new StartEdgeUpdateThread(edges, eIds, this, r.nextInt());
            collection.add(t);
        }

        Pair[] edges = new Pair[numEdgesPerThread + remain];
        int[] eIds = new int[numEdgesPerThread + remain];
        for (int j = 0; j < edges.length; j++) {
            edges[j] = g.getEndpoints(counter);
            eIds[j] = counter;
            counter++;
        }
        EdgeRefinementThread t = new EdgeRefinementThread(edges, eIds, this, numUpdatesPerT, r.nextInt());
        //StartEdgeUpdateThread t = new StartEdgeUpdateThread(edges, eIds, this, r.nextInt());

        collection.add(t);

        try {
            executor.invokeAll(collection);
            //System.out.println("finished");
        } catch (Exception e) {
            e.printStackTrace();
        }

        executor.shutdown();

        double endTime = System.currentTimeMillis();

        gk.setCoord(coord);
        double aktCost = gk.mdlFunctionSimpleSigmoid(); //sigma is calculated but set to 1 again in reInit()
//        if (verbose) {
//            System.out.println(t + " " + aktCost);
//        }
        if (aktCost < bestCost) {
            bestCoord = coord.clone();
            bestCost = aktCost;
        }

        if (display) {
            Visualization v = new Visualization(g);
            v.displayCoordNew(coord, "refine");
            attributedEmbedding.DataObject[] ds = new attributedEmbedding.DataObject[n];
            for (int i = 0; i < n; i++) {
                ds[i] = new attributedEmbedding.DataObject(coord[i], clusterID[i]);
            }
            VisuAttributes vi = new VisuAttributes(ds, "bla");
            vi.setSize(600, 600);
            vi.setLocation(200, 0);
            vi.setVisible(true);
        }
        if (computeCost) {
            System.out.println("cost after refinement: " + bestCost + " " + gk.s.sigma);

        }
        if (verbose) {
            NumberFormat formatter = new DecimalFormat("#0.00000");
            System.out.println("runtime " + formatter.format((endTime - startTime) / 1000d) + " seconds");

        }
        if (writeResult) {
            DataUtils du = new DataUtils();
            du.saveAsMatlab(coord, "refine", "refine.mat");
        }

    }

    public void getBestInitalization(int numTry, int numUpdatesPerEdge, int numT) {
        double startTime = System.currentTimeMillis();
        int numUpdatesPerT = (g.getEdgeCount() / numT) * numUpdatesPerEdge;
        if (verbose) {
            System.out.println(numUpdatesPerT + " updates per thread.");
        }
        int t = 0;
        firstRound(numT);
        Random localRand = new Random(29);
        runUpdates(numT, numUpdatesPerT);
        bestCoord = coord.clone();
        if (computeCost) {
            gk.setCoord(coord);
            bestCost = gk.mdlFunctionSimpleSigmoid();
        } else {
            bestCost = -Double.MAX_VALUE;
        }
        if (verbose) {
            System.out.println(t + " " + bestCost);
        }
        while (t < numTry) {
            t++;
            reInit(localRand.nextInt());
            firstRound(numT);
            runUpdates(numT, numUpdatesPerT);
            gk.setCoord(coord);
            double aktCost = gk.mdlFunctionSimpleSigmoid(); //sigma is calculated but set to 1 again in reInit()
            if (verbose) {
                System.out.println(t + " " + aktCost);
            }
            if (aktCost < bestCost) {
                bestCoord = coord.clone();
                bestCost = aktCost;
            }
        }
        double endTime = System.currentTimeMillis();

        if (display) {
            Visualization v = new Visualization(g);
            v.displayCoordNew(bestCoord, "reInit");
        }
        if (computeCost) {
            System.out.println("cost after init: " + bestCost + " " + gk.s.sigma);

        }
        if (verbose) {
            NumberFormat formatter = new DecimalFormat("#0.00000");
            System.out.println("runtime " + formatter.format((endTime - startTime) / 1000d) + " seconds");

        }
        if (writeResult) {
            DataUtils du = new DataUtils();
            du.saveAsMatlab(bestCoord, "init", "init.mat");
        }
    }

    public void run(int numT, int iter) {
        firstRound(numT);
        //boolean allFree = checkLocks();
        Visualization v = new Visualization(g);
        for (int i = 0; i < iter; i++) {
            runUpdates(numT, 100);
            // boolean allFree = checkLocks(); 
            if (verbose) {
                if (i % 100 == 0) {
                    DataUtils du = new DataUtils();
                    //du.saveAsMatlab(coord, "coord", "coord_"+iter+".mat");
                    double[][] denomV = new double[n][1];
                    for (int j = 0; j < n; j++) {
                        denomV[j][0] = denomVertex[j];
                    }

                    String s = new Integer(i).toString();
                    //du.saveAsMatlab3(coord, nomVertex, denomV, "coord", "nomVertex", "denomVertex", s + ".mat");
                    v.displayCoordNew(coord, s);
                    gk.setCoord(coord);
                    System.out.println(i + " " + gk.mdlFunctionSimpleSigmoid() + " " + gk.s.sigma);
                }
            }
        }
        v.displayCoordNew(coord, "finished");
        gk.setCoord(coord);
        bestCost = gk.mdlFunctionSimpleSigmoid();
        System.out.println("cost: " + bestCost + " " + gk.s.sigma);

    }

    public void runUpdates(int numT, int iterPerT) {
        //public EdgeUpdateThread(Pair<Integer>[] edges, int[] eIds, IncrementalParallel p, int maxIter, int seed) {
        int numEdgesPerThread = g.getEdgeCount() / numT;
        int remain = g.getEdgeCount() - (numT * numEdgesPerThread);
        int counter = 0;
        ExecutorService executor = Executors.newFixedThreadPool(numT);
        Collection collection = new ArrayList();

        for (int i = 0; i < numT - 1; i++) {
            Pair[] edges = new Pair[numEdgesPerThread];
            int[] eIds = new int[numEdgesPerThread];
            for (int j = 0; j < edges.length; j++) {
                edges[j] = g.getEndpoints(counter);
                eIds[j] = counter;
                counter++;
            }
            EdgeUpdateThread t = new EdgeUpdateThread(edges, eIds, this, iterPerT, r.nextInt());
            //StartEdgeUpdateThread t = new StartEdgeUpdateThread(edges, eIds, this, r.nextInt());
            collection.add(t);
        }

        Pair[] edges = new Pair[numEdgesPerThread + remain];
        int[] eIds = new int[numEdgesPerThread + remain];
        for (int j = 0; j < edges.length; j++) {
            edges[j] = g.getEndpoints(counter);
            eIds[j] = counter;
            counter++;
        }
        EdgeUpdateThread t = new EdgeUpdateThread(edges, eIds, this, iterPerT, r.nextInt());
        //StartEdgeUpdateThread t = new StartEdgeUpdateThread(edges, eIds, this, r.nextInt());

        collection.add(t);

        try {
            executor.invokeAll(collection);
            //System.out.println("finished");
        } catch (Exception e) {
            e.printStackTrace();
        }

        executor.shutdown();

    }

    //numTry: number of attempts to sample a distant object, similar as in FastMap
    private void initClusters(int numClusters, int numTry, int numIter) {
        this.k = numClusters;
        this.clusters = new Cluster[k];
        for (int i = 0; i < k; i++) {
            clusters[i] = new Cluster(d);
        }
        clusterID = new int[n];

        Random cl = new Random(22);
        //sample first centroid randomly
        int c1 = r.nextInt(n);
        boolean[] selected = new boolean[n];
        selected[c1] = true;
        int count = 1;
        Vector<Integer> indexSelected = new Vector<Integer>();
        indexSelected.add(c1);
        //select distant centroids
        while (count < numClusters) {
            double[] sumDist = new double[numTry];
            int[] index = new int[numTry];
            double bestSumDist = 0.0;
            int bestIndex = -1;
            for (int i = 0; i < numTry; i++) {
                int cNext = r.nextInt(n);
                while (selected[cNext]) {
                    cNext = r.nextInt(n);
                }
                selected[cNext] = true;
                index[i] = cNext;
                for (int j = 0; j < indexSelected.size(); j++) {
                    sumDist[i] += dist(coord[cNext], coord[indexSelected.elementAt(j)]);
                }
            }
            for (int i = 0; i < numTry; i++) {
                if (sumDist[i] > bestSumDist) {
                    bestSumDist = sumDist[i];
                    bestIndex = index[i];
                }
            }
            indexSelected.add(bestIndex);
            selected[bestIndex] = true;
            count++;
        }
        //perform first assignment
        Object[] bla = indexSelected.toArray();
        Integer[] initialCenters = new Integer[bla.length];
        for (int i = 0; i < bla.length; i++) {
            initialCenters[i] = (Integer) bla[i];
        }
        for (int i = 0; i < n; i++) {
            double minDist = Double.MAX_VALUE;
            int minIndex = -1;
            for (int j = 0; j < k; j++) {
                double aktDist = dist(coord[initialCenters[j]], coord[i]);
                if (aktDist < minDist) {
                    minDist = aktDist;
                    minIndex = j;
                }
            }
            clusterID[i] = minIndex;
            clusters[minIndex].members.add(i);
            for (int j = 0; j < d; j++) {
                clusters[minIndex].sum[j] += coord[i][j];
            }

        }

        if(display){
        attributedEmbedding.DataObject[] ds = new attributedEmbedding.DataObject[n];
        for (int i = 0; i < n; i++) {
            ds[i] = new attributedEmbedding.DataObject(coord[i], clusterID[i]);
        }
        VisuAttributes vi = new VisuAttributes(ds, "first assignment");
        vi.setSize(600, 600);
        vi.setLocation(200, 0);
        vi.setVisible(true);
        }
        runKmeans(numIter);

    }

    private void runKmeans(int numIter) {
        for (int i = 0; i < numIter; i++) {
            kmeansIteration();
        }
    }

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }

    private void resetClusters() {
        for (int i = 0; i < k; i++) {
            clusters[i] = new Cluster(d);
        }
    }

    //round-wise update and assignment
    private void kmeansIteration() {
        double[][] means = new double[k][d];
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < d; j++) {
                means[i][j] = clusters[i].sum[j] / clusters[i].members.size();
            }
        }
        resetClusters();
        for (int i = 0; i < n; i++) {
            double minDist = Double.MAX_VALUE;
            int minIndex = -1;
            for (int j = 0; j < k; j++) {
                double aktDist = dist(means[j], coord[i]);
                if (aktDist < minDist) {
                    minDist = aktDist;
                    minIndex = j;
                }
            }
            clusterID[i] = minIndex;
            clusters[minIndex].members.add(i);
            for (int j = 0; j < d; j++) {
                clusters[minIndex].sum[j] += coord[i][j];
            }
        }

    }

    private void update() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
