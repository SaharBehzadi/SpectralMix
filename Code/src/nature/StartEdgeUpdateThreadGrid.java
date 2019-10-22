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
public class StartEdgeUpdateThreadGrid implements Callable {

    Pair<Integer>[] edges; //the subset of edges this thread updates
    int[] eIds; //ids of the edges
    IncrementalParallelGrid p;
    Random r;
    int n;
    int d;
    int iter;
    int maxIter;

    public StartEdgeUpdateThreadGrid(Pair<Integer>[] edges, int[] eIds, IncrementalParallelGrid p, int seed) {
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
        p.locks[Math.min(edge.getFirst(), edge.getSecond())].lock();
        p.locks[Math.max(edge.getFirst(), edge.getSecond())].lock();
//perform the update
        try {
            double x = dist(p.coord[edge.getFirst()], p.coord[edge.getSecond()]);
            double[] dw = p.s.parabola(x, p.isEdge[i]);
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
//end of update
        } finally {
//release the locks
            p.locks[edge.getFirst()].unlock();
            p.locks[edge.getSecond()].unlock();
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
