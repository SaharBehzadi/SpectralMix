/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.util.Pair;
import java.util.Random;

/**
 *
 * @author plantc59cs
 */
public class GempeParNode implements Runnable {

    int iterations; //number of updates this thread will perform
    Pair<Integer>[] edges; //the subset of edges this thread updates
    GempePar p;
    int n;
    int d;
    Random r;
    int seed;
    boolean verbose = true;

    public GempeParNode(int iterations, Pair<Integer>[] edges, GempePar p, int seed) {
        this.iterations = iterations;
        this.edges = edges;
        this.p = p;
        this.seed = seed;
        n = p.g.getVertexCount();
        d = p.coord[0].length;
        r = new Random(seed);

    }

    public void run() {
        if (verbose) {
            System.out.println("thread with seed " + seed + " started. Proessing " + edges.length + " edges.");
        }
        for (int i = 0; i < iterations; i++) {
            double[] sumWeights = new double[n];
            double[][] positions = new double[n][d];
            for (int j = 0; j < edges.length; j++) {
                //update the positions of both endpoints of this edge
                
                
                //first the edge itself
                double x = dist(p.coord[edges[j].getFirst()], p.coord[edges[j].getSecond()]);
                double[] dw = p.s.parabola(x, true);
                sumWeights[edges[j].getFirst()] += dw[1];
                sumWeights[edges[j].getSecond()] += dw[1];
                double sij = 0;
                for (int k = 0; k < d; k++) {
                    sij += Math.pow((p.coord[edges[j].getFirst()][k] - p.coord[edges[j].getSecond()][k]), 2);
                }
                if (sij != 0) {
                    sij = dw[0] / Math.sqrt(sij);
                }

                for (int k = 0; k < d; k++) {
                    positions[edges[j].getFirst()][k] += dw[1] * (p.coord[edges[j].getSecond()][k] + sij * (p.coord[edges[j].getFirst()][k] - p.coord[edges[j].getSecond()][k]));
                    positions[edges[j].getSecond()][k] += dw[1] * (p.coord[edges[j].getFirst()][k] + sij * (p.coord[edges[j].getSecond()][k] - p.coord[edges[j].getFirst()][k]));
                }

                int notEdge1 = r.nextInt(n);
                while (p.g.isNeighbor(edges[j].getFirst(), notEdge1)) {
                    notEdge1 = r.nextInt(n);
                }
                x = dist(p.coord[edges[j].getFirst()], p.coord[notEdge1]);
                dw = p.s.parabola(x, false);
                sumWeights[edges[j].getFirst()] += dw[1];
                sumWeights[notEdge1] += dw[1];
                sij = 0;
                for (int k = 0; k < d; k++) {
                    sij += Math.pow((p.coord[edges[j].getFirst()][k] - p.coord[notEdge1][k]), 2);
                }
                if (sij != 0) {
                    sij = dw[0] / Math.sqrt(sij);
                }
                for (int k = 0; k < d; k++) {
                    positions[edges[j].getFirst()][k] += dw[1] * (p.coord[notEdge1][k] + sij * (p.coord[edges[j].getFirst()][k] - p.coord[notEdge1][k]));
                    positions[notEdge1][k] += dw[1] * (p.coord[edges[j].getFirst()][k] + sij * (p.coord[notEdge1][k] - p.coord[edges[j].getFirst()][k]));
                }

                int notEdge2 = r.nextInt(n);
                while (p.g.isNeighbor(edges[j].getSecond(), notEdge2)) {
                    notEdge2 = r.nextInt(n);
                }
                x = dist(p.coord[edges[j].getSecond()], p.coord[notEdge2]);
                dw = p.s.parabola(x, false);
                sumWeights[edges[j].getSecond()] += dw[1];
                sumWeights[notEdge2] += dw[1];
                sij = 0;
                for (int k = 0; k < d; k++) {
                    sij += Math.pow((p.coord[edges[j].getSecond()][k] - p.coord[notEdge2][k]), 2);
                }
                if (sij != 0) {
                    sij = dw[0] / Math.sqrt(sij);
                }
                for (int k = 0; k < d; k++) {
                    positions[edges[j].getSecond()][k] += dw[1] * (p.coord[notEdge2][k] + sij * (p.coord[edges[j].getSecond()][k] - p.coord[notEdge2][k]));
                    positions[notEdge2][k] += dw[1] * (p.coord[edges[j].getSecond()][k] + sij * (p.coord[notEdge2][k] - p.coord[edges[j].getSecond()][k]));
                }
            }//edges
            for (int l = 0; l < n; l++) {
                for (int k = 0; k < d; k++) {
                    //sumWeights[i] = Math.max(sumWeights[i], 1E-9);
                    if (!(positions[l][k] == 0 && sumWeights[l] == 0)) {
                        //   System.out.println("m");
                        p.coord[l][k] = positions[l][k] / sumWeights[l];
                    }
                }
            }

        }
        if (verbose) {
            System.out.println(iterations + " iterations finished.");
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

}
