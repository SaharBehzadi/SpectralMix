/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.util.Pair;
import java.util.Random;
import java.util.concurrent.Callable;

/**
 *
 * @author plantc59cs
 */
public class StartEdgeUpdateThread implements Callable {

    Pair<Integer>[] edges; //the subset of edges this thread updates
    int[] eIds; //ids of the edges
    IncrementalParallel p;
    Random r;
    int n;
    int d;
    int iter;
    int maxIter;

    public StartEdgeUpdateThread(Pair<Integer>[] edges, int[] eIds, IncrementalParallel p, int seed) {
        this.edges = edges;
        this.eIds = eIds;
        this.p = p;
        this.maxIter = maxIter;
        this.n = p.n;
        this.d = p.d;
        iter = 0;
        r = new Random(seed);
    }

    public String call() {
        for (int i = 0; i < edges.length; i++) {
            //System.out.println(i);
            update(edges[i], eIds[i]);
        }
        return ("finished");

    }

    private void update(Pair<Integer> edge, int i) {
        //get locks for the two endpoints
//        p.locks[Math.min(edge.getFirst(), edge.getSecond())].lock();
//        p.locks[Math.max(edge.getFirst(), edge.getSecond())].lock();
        p.lock4u(edge.getFirst(), edge.getSecond(), p.vNonEdge1[i], p.vNonEdge2[i]);
        int formerNonEdge1 = p.vNonEdge1[i];
        int formerNonEdge2 = p.vNonEdge2[i];
        //try to get locks for not connected endpoints
        int notEdge1 = getNonEdge1(edge);
        boolean free1 = p.locks[notEdge1].tryLock();
        while (!free1) {
            notEdge1 = getNonEdge1(edge);
            free1 = p.locks[notEdge1].tryLock();
        }
        int notEdge2 = getNonEdge2(edge, notEdge1);
        boolean free2 = p.locks[notEdge2].tryLock();
        while (!free2) {
            notEdge2 = getNonEdge2(edge, notEdge1);
            free2 = p.locks[notEdge2].tryLock();
        }
//        boolean freeFirst = p.locks[edge.getFirst()].isLocked();
//        boolean freeSecond = p.locks[edge.getSecond()].isLocked();
//        boolean freeNotEdge1 = p.locks[notEdge1].isLocked();
//        boolean freeNotEdge2 = p.locks[notEdge2].isLocked();
//perform the update
        try {

//            //            //DEBUG
//            double[] dwT = p.s.parabola(2.0, false);
//            double[] dwTT = p.s.parabola(5.0, false);
//            //DEBUG
            double x = dist(p.coord[edge.getFirst()], p.coord[edge.getSecond()]);
            double[] dw = p.s.parabola(x, true);
            double[] h = zij(dw[0], dw[1], edge.getFirst(), edge.getSecond());
            for (int j = 0; j < d; j++) {
                p.nomVertex[edge.getFirst()][j] += h[j] - (p.valid[i] ? p.nomEdge[i][j] : 0);
            }
            p.nomEdge[i] = h;
            p.denomVertex[edge.getFirst()] += dw[1] - (p.valid[i] ? p.denomEdge[i] : 0);
            h = zij(dw[0], dw[1], edge.getSecond(), edge.getFirst());
            for (int j = 0; j < d; j++) {
                p.nomVertex[edge.getSecond()][j] += h[j] - (p.valid[i] ? p.nomEdgeR[i][j] : 0);
            }
            p.nomEdgeR[i] = h;
            p.denomVertex[edge.getSecond()] += dw[1] - (p.valid[i] ? p.denomEdge[i] : 0);
            p.denomEdge[i] = dw[1];

            x = dist(p.coord[edge.getFirst()], p.coord[notEdge1]);
            dw = p.s.parabola(x, false);

//            //DEBUG
//            double[] dwT = p.s.parabola(0.2, true);
//            double[] dwTT = p.s.parabola(0.3, true);
//            //DEBUG
            h = zij(dw[0], dw[1], edge.getFirst(), notEdge1);
            for (int j = 0; j < d; j++) {
                p.nomVertex[edge.getFirst()][j] += h[j] - (p.valid[i] ? p.nomNonEdge1[i][j] : 0);
            }
            p.nomNonEdge1[i] = h;
            p.denomVertex[edge.getFirst()] += dw[1] - (p.valid[i] ? p.denomNonEdge1[i] : 0);
            h = zij(dw[0], dw[1], notEdge1, edge.getFirst());
            for (int j = 0; j < d; j++) {
                p.nomVertex[notEdge1][j] += h[j];// - (valid[i] ? nomNonEdge1R[i][j] : 0);
                if (p.valid[i]) {
                    p.nomVertex[p.vNonEdge1[i]][j] -= p.nomNonEdge1R[i][j];
                }
            }
            p.nomNonEdge1R[i] = h;
            p.denomVertex[notEdge1] += dw[1];// - (valid[i] ? denomNonEdge1[i] : 0);
            if (p.valid[i]) {
                p.denomVertex[p.vNonEdge1[i]] -= p.denomNonEdge1[i];
            }
            p.denomNonEdge1[i] = dw[1];
            p.vNonEdge1[i] = notEdge1;

            x = dist(p.coord[edge.getSecond()], p.coord[notEdge2]);
            //double aktW = ((2 * s.mu - 2 * x) * (s.a - s.b)) / (Math.pow(Math.PI, 0.5) * Math.pow(s.sigma, 3) * Math.exp(Math.pow(s.mu - x, 2) / Math.pow(s.sigma, 2)));
            dw = p.s.parabola(x, false);
            h = zij(dw[0], dw[1], edge.getSecond(), notEdge2);
            for (int j = 0; j < d; j++) {
                p.nomVertex[edge.getSecond()][j] += h[j] - (p.valid[i] ? p.nomNonEdge2[i][j] : 0);
            }
            p.nomNonEdge2[i] = h;
            p.denomVertex[edge.getSecond()] += dw[1] - (p.valid[i] ? p.denomNonEdge2[i] : 0);
            h = zij(dw[0], dw[1], notEdge2, edge.getSecond());
            for (int j = 0; j < d; j++) {
                p.nomVertex[notEdge2][j] += h[j];
                if (p.valid[i]) {
                    p.nomVertex[p.vNonEdge2[i]][j] -= p.nomNonEdge2R[i][j];
                }
            }
            p.nomNonEdge2R[i] = h;
            p.denomVertex[notEdge2] += dw[1];
            if (p.valid[i]) {
                p.denomVertex[p.vNonEdge2[i]] -= p.denomNonEdge2[i];
            }
            p.denomNonEdge2[i] = dw[1];
            p.vNonEdge2[i] = notEdge2;

            p.valid[i] = true;

//end of update
        } finally {
//release the locks
//            freeFirst = p.locks[edge.getFirst()].isLocked();
//            freeSecond = p.locks[edge.getSecond()].isLocked();
//            freeNotEdge1 = p.locks[notEdge1].isLocked();
//            freeNotEdge2 = p.locks[notEdge2].isLocked();

            p.locks[edge.getFirst()].unlock();
            p.locks[edge.getSecond()].unlock();
            p.locks[formerNonEdge1].unlock();
            p.locks[formerNonEdge2].unlock();
            p.locks[notEdge1].unlock();
            p.locks[notEdge2].unlock();

//            freeFirst = p.locks[edge.getFirst()].isLocked();
//            freeSecond = p.locks[edge.getSecond()].isLocked();
//            freeNotEdge1 = p.locks[notEdge1].isLocked();
//            freeNotEdge2 = p.locks[notEdge2].isLocked();
        }
    }

    private double[] zij(double dij, double wij, int i, int j) {
        double[] res = new double[d];
        double sij = 0.0;
        for (int k = 0; k < d; k++) {
            sij += Math.pow((p.coord[i][k] - p.coord[j][k]), 2);
        }
        if (sij != 0) {
            sij = dij / Math.sqrt(sij);
        }
        for (int k = 0; k < d; k++) {
            res[k] = wij * (p.coord[j][k] + sij * (p.coord[i][k] - p.coord[j][k]));
        }
        return res;
    }

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }

    private int getNonEdge1(Pair<Integer> edge) {
        int nonEdge1 = r.nextInt(n);
        while (p.g.isNeighbor(edge.getFirst(), nonEdge1)) {
            nonEdge1 = r.nextInt(n);
        }
        return nonEdge1;
    }

    private int getNonEdge2(Pair<Integer> edge, int nonEdge1) {
        int nonEdge2 = r.nextInt(n);
        while (p.g.isNeighbor(edge.getSecond(), nonEdge2) || nonEdge1 == nonEdge2) {
            nonEdge2 = r.nextInt(n);
        }
        return nonEdge2;
    }

}
