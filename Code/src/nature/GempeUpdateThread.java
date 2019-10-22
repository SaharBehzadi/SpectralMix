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
public class GempeUpdateThread implements Runnable {

    int iterations; //number of updates this thread will perform
    Pair<Integer>[] edges; //the subset of edges this thread updates
    GempePar p;
    int n;
    int d;
    Random r;
    DataUtils du;
    int seed;
    boolean verbose = false;

    public GempeUpdateThread(int iterations, Pair<Integer>[] edges, GempePar p, int seed) {
        this.iterations = iterations;
        this.edges = edges;
        this.p = p;
        this.seed = seed;
        n = p.g.getVertexCount();
        d = p.coord[0].length;
        r = new Random(seed);
        du = new DataUtils();

    }

    public void run() {
        if (verbose) {
            System.out.println("thread with seed " + seed + " started. Processing " + edges.length + " edges.");
        }
        for (int i = 0; i < iterations; i++) {
           // r = new Random(seed + i);
            for (int j = 0; j < edges.length; j++) {
//                if (j == 3 && i == 7) {
//                    System.out.println(j);
//                }
                update(edges[j], true);
                int notEdge1 = r.nextInt(n);
                
                while (p.g.isNeighbor(edges[j].getFirst(), notEdge1) || edges[j].getFirst() == notEdge1) {
                    notEdge1 = r.nextInt(n);
                }
                Pair<Integer> e = new Pair<Integer>(edges[j].getFirst(), notEdge1);
                update(e, false);
                int notEdge2 = r.nextInt(n);
                
                while (p.g.isNeighbor(edges[j].getSecond(), notEdge2) || edges[j].getSecond() == notEdge2) {
                    notEdge2 = r.nextInt(n);
                }
                e = new Pair<Integer>(edges[j].getSecond(), notEdge2);
                update(e, false);
            }
        }
        if (verbose) {
            System.out.println(iterations + " iterations finished.");
        }
//        p.computeSigmoid();
//        Visualization v = new Visualization(p.g);
//        v.displayCoordNew(p.coord, Integer.toString(seed));

    }

    public void update(Pair<Integer> edge, boolean connected) {
        double x = dist(p.coord[edge.getFirst()], p.coord[edge.getSecond()]);
        double[] dw = p.s.parabola(x, connected);
        double[] coordNewFirst = new double[d];
        double[] coordNewSecond = new double[d];
        int index_ij = du.getIndex(edge.getFirst(), edge.getSecond(), n);
        double nodeWeightNewFirst = p.nodeWeights[edge.getFirst()] + dw[0] - p.edgeWeights[index_ij];
        double nodeWeightNewSecond = p.nodeWeights[edge.getSecond()] + dw[0] - p.edgeWeights[index_ij];
        double sij = 0;
        for (int k = 0; k < d; k++) {
            sij += Math.pow((p.coord[edge.getFirst()][k] - p.coord[edge.getSecond()][k]), 2);
        }
        if (sij != 0) {
            sij = dw[0] / Math.sqrt(sij);
        }
        double[] zijNew = new double[d];
        double[] zjiNew = new double[d];
        for (int k = 0; k < d; k++) {
            zijNew[k] = dw[1] * (p.coord[edge.getSecond()][k] + sij * (p.coord[edge.getFirst()][k] - p.coord[edge.getSecond()][k]));
            coordNewFirst[k] = (p.nodeWeights[edge.getFirst()] * p.coord[edge.getFirst()][k] - p.z[edge.getFirst()][edge.getSecond()][k] + zijNew[k]) / nodeWeightNewFirst;
        }
        for (int k = 0; k < d; k++) {
            zjiNew[k] = dw[1] * (p.coord[edge.getFirst()][k] + sij * (p.coord[edge.getSecond()][k] - p.coord[edge.getFirst()][k]));
            coordNewSecond[k] = (p.nodeWeights[edge.getSecond()] * p.coord[edge.getSecond()][k] - p.z[edge.getFirst()][edge.getSecond()][k] + zjiNew[k]) / nodeWeightNewSecond;
        }
        p.edgeWeights[index_ij] = dw[1];
        p.nodeWeights[edge.getFirst()] = nodeWeightNewFirst;
        p.nodeWeights[edge.getSecond()] = nodeWeightNewSecond;
        p.coord[edge.getFirst()] = coordNewFirst;
        p.coord[edge.getSecond()] = coordNewSecond;
        p.z[edge.getFirst()][edge.getSecond()] = zijNew;
        p.z[edge.getSecond()][edge.getFirst()] = zjiNew;

    }

    public void runOld() {
        if (verbose) {
            System.out.println("thread with seed " + seed + " started. Proessing " + edges.length + " edges.");
        }
        DataUtils du = new DataUtils();
        for (int i = 0; i < iterations; i++) {
            double[] sumWeights = new double[n];
            double[][] positions = new double[n][d];
            for (int j = 0; j < edges.length; j++) {
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
