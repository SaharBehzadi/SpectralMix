/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import java.awt.geom.Point2D;
import java.util.List;
import java.util.Random;
import mdsj.MDSJ;
import mdsj.StressMinimization;

/**
 *
 * @author claudia
 */
public class Embedding {

    Graph g;
    int dim;
    double[][] coord;
    int[] clid;
    int numObj;

    public Embedding(Graph g, int dim) {
        this.g = g;
        this.dim = dim;
        numObj = g.getVertexCount();
        coord = new double[numObj][dim];

    }

    public void setCoord(double[][] coord) {
        this.coord = coord;
        GraphCompression gk = new GraphCompression(g, coord, 30);
        double mdlHist = gk.mdlHisto();
        Visualization v = new Visualization(g);
        v.displayCoord(coord, "init: " + mdlHist);

    }

    public void setClid(int[] clid) {
        this.clid = clid;
    }

    public void randomInit(int randomSeed) {
        Random r = new Random(randomSeed);
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < coord[i].length; j++) {
                coord[i][j] = r.nextDouble();
            }
        }
        GraphCompression gk = new GraphCompression(g, coord, 30);
        double mdlHist = gk.mdlHisto();
        Visualization v = new Visualization(g);
        v.displayCoordNew(coord, "random init: " + mdlHist, clid);
        //      v.displayCoord(coord, "random init: " );
    }

    public void distanceScaling(int iter) {
        double[][] dist = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j) || (i == j)) {
                    dist[i][j] = 0.0;
                } else {
                    dist[i][j] = 1.0;
                }
            }
        }
        double[][] weights = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                weights[i][j] = 1.0;
            }
        }
        double[][] c = new Matrix(coord).transpose().getArrayCopy();
        StressMinimization st = new StressMinimization(dist, c, weights);
        StressMinimization.majorize(c, dist, weights, iter);
//            System.out.println(st.iterate(5000));
        double[][] cc = st.getPositions();
        coord = new Matrix(cc).transpose().getArrayCopy();
        GraphCompression gkk = new GraphCompression(g, coord, 30);
        double mdlHistt = gkk.mdlHisto();
        Visualization v = new Visualization(g);
        v.displayCoord(coord, "majorize " + mdlHistt);

    }

    public void distanceScalingUpdateWeights(int iter) {
        double[][] dist = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j) || (i == j)) {
                    dist[i][j] = 0.0;
                } else {
                    dist[i][j] = 1.0;
                }
            }
        }
        double[][] weights = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                weights[i][j] = 1.0;
            }
        }
        double[][] c = new Matrix(coord).transpose().getArrayCopy();
        StressMinimization st = new StressMinimization(dist, c, weights);
        StressMinimization.majorize(c, dist, weights, 1);
        int ii = 0;
        while (ii < iter) {
        }


//            System.out.println(st.iterate(5000));
        double[][] cc = st.getPositions();
        coord = new Matrix(cc).transpose().getArrayCopy();
        GraphCompression gkk = new GraphCompression(g, coord, 30);
        double mdlHistt = gkk.mdlHisto();
        Visualization v = new Visualization(g);
        v.displayCoord(coord, "majorize " + mdlHistt);

    }

    public void mdsInit() {
        double[][] dist = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j) || (i == j)) {
                    dist[i][j] = 0.0;
                } else {
                    dist[i][j] = 1.0;
                }
            }
        }

        Visualization v = new Visualization(g);
        double[][] init = MDSJ.classicalScaling(dist); // apply MDS
        coord = new Matrix(init).transpose().getArrayCopy();
        GraphCompression gk = new GraphCompression(g, coord, 30);
        double mdlHist = gk.mdlHisto();
        v.displayCoord(coord, "mds: " + mdlHist);

    }

    public Point2D[] getLaplacianCoord() {
        Point2D.Double[] coord = new Point2D.Double[numObj];
        double[][] w = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j)) {
                    w[i][j] = 1.0;
                } else {
                    w[i][j] = 0.0;
                }
            }
        }
        double[][] d = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                d[i][i] += w[i][j];
            }
        }

//        //laplacian eigenmaps
        Matrix D = new Matrix(d);
        Matrix LL = D.minus(new Matrix(w));
        //debug generalized eigenvalue problem: write out LL and D
        IO ea = new IO();
        ea.writeDoubleToMatlab(d, "D");
        ea.writeDoubleToMatlab(LL.getArrayCopy(), "L");

        double[][] dd = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            dd[i][i] = 1.0 / d[i][i];
        }
        Matrix DD = new Matrix(dd);
        Matrix Ln = DD.times(LL);


        EigenvalueDecomposition eig = new EigenvalueDecomposition(new Matrix(w));
        Matrix vectors = eig.getV();
        Matrix values = eig.getD();
        //debug
        ea.writeDoubleToMatlab(vectors.getArrayCopy(), "vectors");
        ea.writeDoubleToMatlab(values.getArrayCopy(), "values");

        for (int i = 0; i < numObj; i++) {
            coord[i] = new Point2D.Double();
            //2 largest eigenvalues
//            double x = vectors.get((vectors.getRowDimension()-1), i) * values.get(values.getRowDimension()-1, values.getRowDimension()-1);
//            double y = vectors.get((vectors.getRowDimension()-2), i) * values.get(values.getRowDimension()-2, values.getRowDimension()-2);

            //2 smallest eigenvalues > 0

            coord[i].setLocation(vectors.get(1, i), vectors.get(2, i));
        }
        return coord;

    }




    public Point2D[] getIsoCoord() {
        Point2D[] pp = new Point2D[g.getVertexCount()];

        DijkstraShortestPath<Integer, Integer> alg = new DijkstraShortestPath(g);
        double[][] distn = new double[numObj][numObj];
        double maxDist_mds = 0.0;
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (i > j) {
                    List<Integer> l = alg.getPath(i, j);
                    if (l.size() > 0) {
                        distn[i][j] = l.size();
                        distn[j][i] = distn[i][j];
                        if (distn[i][j] > maxDist_mds) {
                            maxDist_mds = distn[i][j];
                        }
                    }
                }
            }
        }
        for (int i = 0; i < distn.length; i++) {
            for (int j = 0; j < distn.length; j++) {
                if (distn[i][j] == 0.0 && i != j) {
                    distn[i][j] = maxDist_mds;
                }
            }
        }

        double[][] init = MDSJ.classicalScaling(distn); // apply MDS

        double[][] c = new Matrix(init).transpose().getArrayCopy();
        for (int i = 0; i < coord.length; i++) {
            pp[i] = new Point2D.Double();
            pp[i].setLocation(c[i][0], c[i][1]);
        }
        return pp;
    }

    public void isoInit() {

        DijkstraShortestPath<Integer, Integer> alg = new DijkstraShortestPath(g);
        double[][] distn = new double[numObj][numObj];
        double maxDist_mds = 0.0;
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (i > j) {
                    List<Integer> l = alg.getPath(i, j);
                    if (l.size() > 0) {
                        distn[i][j] = l.size();
                        distn[j][i] = distn[i][j];
                        if (distn[i][j] > maxDist_mds) {
                            maxDist_mds = distn[i][j];
                        }
                    }
                }
            }
        }
        for (int i = 0; i < distn.length; i++) {
            for (int j = 0; j < distn.length; j++) {
                if (distn[i][j] == 0.0 && i != j) {
                    distn[i][j] = maxDist_mds;
                }
            }
        }
        Visualization v = new Visualization(g);
        double[][] init = MDSJ.classicalScaling(distn); // apply MDS
        GraphCompression gk = new GraphCompression(g, new Matrix(init).transpose().getArrayCopy(), 30);
        double mdlHist = gk.mdlHisto();
        double[][] c = new Matrix(init).transpose().getArrayCopy();
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < dim; j++) {
                coord[i][j] = c[i][j];
            }
        }
        //   System.arraycopy(new Matrix(init).transpose().getArrayCopy(), 0, coord, 0, dim);
        v.displayCoordNew(coord, "iso: " + mdlHist, clid);
        //v.displayCoord(coord, "iso: ");

    }

    public void ownInit() {

        double[][] bestCoord = new double[numObj][dim];
        int numVertices = g.getVertexCount();
        int numEdges = g.getEdgeCount();
        int numPossEdges = (numVertices * (numVertices - 1)) / 2;
        double pLink_g = (double) (numEdges * 2) / (double) numVertices;
        double pNoLink_g = 1.0 - pLink_g;
        GraphCompression gc = new GraphCompression(g, coord, 30);
        System.out.println("codingCostNoEmbedding: " + gc.codingCostEmbedding());
        double baseMDL = 0.0;
        baseMDL = gc.mdlHisto();
        System.out.println("baseMDL: " + baseMDL);
        int iter = 0;
        double bestMdl = java.lang.Double.MAX_VALUE;
        while (iter < 100) {
            for (int i = 0; i < numObj; i++) {
                // System.out.println(i);
                double[] meanLink = new double[dim];
                double[] meanNotLink = new double[dim];
                int counter_links = 0;
                int counter_noLinks = 0;
                for (int k = 0; k < numObj; k++) {
                    if (g.isNeighbor(i, k)) {
                        counter_links++;
                        for (int j = 0; j < dim; j++) {
                            meanLink[j] += coord[k][j];
                        }
                    } else {
                        counter_noLinks++;
                        for (int j = 0; j < dim; j++) {
                            meanNotLink[j] += coord[k][j];
                        }
                    }
                }
                if (counter_noLinks > 0 && counter_links > 0) {
                    for (int j = 0; j < dim; j++) {
                        meanLink[j] /= (double) counter_links;
                        meanNotLink[j] /= (double) counter_noLinks;
                    }
                    double[] maxDistantPoint = new double[dim];
                    for (int k = 0; k < numObj; k++) {
                        for (int j = 0; j < dim; j++) {
                            if ((Math.abs(coord[k][j] - meanNotLink[j])) > maxDistantPoint[j]) {
                                maxDistantPoint[j] = coord[k][j];
                            }
                        }
                    }

                    if (iter % 2 == 0) {
                        for (int j = 0; j < dim; j++) {
                            coord[i][j] = (coord[i][j] + meanLink[j]) / 2;
                        }
                    } else {

                        //option 2 : move away from not links
                        for (int j = 0; j < dim; j++) {
                            coord[i][j] = (coord[i][j] + maxDistantPoint[j]) / 2;
                        }
                    }

//                        //weighted sum
//                         for (int j = 0; j < dim; j++) {
//                            coord[i][j] = (((double) (counter_noLinks) * maxDistantPoint[j]) + ((double) (counter_links) * coord[i][j])) / (counter_links + counter_noLinks);
//                        }


                }
//                    if (iter % 2 == 0) {
//                        //System.out.println(iter);
//                        System.arraycopy(meanLink, 0, coord[i], 0, dim);
//                    } else {
//                        double[] maxDistantPoint = new double[dim];
//                        for (int k = 0; k < numObj; k++) {
//                            for (int j = 0; j < dim; j++) {
//                                if ((Math.abs(coord[k][j] - meanNotLink[j])) > maxDistantPoint[j]) {
//                                    maxDistantPoint[j] = coord[i][j];
//                                }
//                            }
//                        }
//                        for (int j = 0; j < dim; j++) {
//                            coord[i][j] = ((double) (counter_noLinks) * maxDistantPoint[j]) + ((double) (counter_links) * coord[i][j]) / (counter_links + counter_noLinks);
//                        }
//                        System.arraycopy(maxDistantPoint, 0, coord[i], 0, dim);
//                    }

                if (counter_links == 0 && counter_noLinks > 0) {
                    double[] maxDistantPoint = new double[dim];
                    for (int k = 0; k < numObj; k++) {
                        for (int j = 0; j < dim; j++) {
                            if ((Math.abs(coord[k][j] - meanNotLink[j])) > maxDistantPoint[j]) {
                                maxDistantPoint[j] = coord[k][j];
                            }
                        }
                    }
                    System.arraycopy(maxDistantPoint, 0, coord[i], 0, dim);

                }
                if (counter_links > 0 && counter_noLinks == 0) {
                    for (int j = 0; j < dim; j++) {
                        meanLink[j] /= (double) counter_links;
                    }
                    System.arraycopy(meanLink, 0, coord[i], 0, dim);
                }
            }

            gc = new GraphCompression(g, coord, 30);
//            if (iter == 89) {
//                System.out.println("m");
//            }
//            if  (iter == 95){
//                System.out.println("m");
//            }
            double newMdl = 0.0;

            newMdl = gc.mdlHisto();

            //double newMdl = gc.mdlFunction();
            System.out.println("iter: " + iter + " " + newMdl);
            if (newMdl < bestMdl) {
                bestMdl = newMdl;
                for (int l = 0; l < numObj; l++) {
                    System.arraycopy(coord[l], 0, bestCoord[l], 0, dim);
                }
                //re-scale bestCoordinates between 0 and 1
                double[] min = new double[dim];
                double[] max = new double[dim];
                for (int j = 0; j < dim; j++) {
                    min[j] = java.lang.Double.MAX_VALUE;
                    max[j] = -min[j];
                }
                for (int l = 0; l < numObj; l++) {
                    for (int j = 0; j < dim; j++) {
                        if (bestCoord[l][j] < min[j]) {
                            min[j] = bestCoord[l][j];
                        }
                        if (bestCoord[l][j] > max[j]) {
                            max[j] = bestCoord[l][j];
                        }
                    }
                }
                for (int l = 0; l < numObj; l++) {
                    for (int j = 0; j < dim; j++) {
                        bestCoord[l][j] = (bestCoord[l][j] - min[j]) / (max[j] - min[j]);

                    }
                }

            }
            iter++;

        }
        System.out.println("bestMdl: " + bestMdl);
        Visualization v = new Visualization(g);
        v.displayCoordNew(bestCoord, "own init " + bestMdl, clid);
        System.arraycopy(bestCoord, 0, coord, 0, dim);
        //display best coordinates
    }

    public void improveCoordinates() {
        int n = coord.length;
        int d = coord[0].length;
        int m = n * (n - 1) / 2;
        Random r = new Random(2);
        int count = 0;
        for (double sig = 0.05; sig > 0.00001; sig *= 0.98) {
            if (count % 100 == 0) {
                System.out.println(count);
            }
            count++;
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), g.isNeighbor(i, j));
                }
            }
            Sigmoid s = new Sigmoid(p);
            for (int i = 0; i < n; i++) {
                // now determine individual cost of all edges and non-edges of node i
                double[] newCoord = coord[i].clone();
                double costBefore = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        if (g.isNeighbor(i, j)) {
                            costBefore -= log2(s.f(dist(newCoord, coord[j])));
                        } else {
                            costBefore -= log2(1.0 - s.f(dist(newCoord, coord[j])));
                        }
                    }
                }
                for (int j = 0; j < d; j++) {
                    newCoord[j] += r.nextGaussian() * sig;
                }
                // now determine new cost
                double costAfter = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        if (g.isNeighbor(i, j)) {
                            costAfter -= log2(s.f(dist(newCoord, coord[j])));
                        } else {
                            costAfter -= log2(1.0 - s.f(dist(newCoord, coord[j])));
                        }
                    }
                }
                if (costAfter < costBefore) {
                    coord[i] = newCoord;
                }
            }
            if (sig < 0.00001) {
                s = s;
            }
        }
        GraphCompression gc = new GraphCompression(g, coord, 30);
        double newMdl = gc.mdlHisto();
        Visualization v = new Visualization(g);
        v.displayCoordNew(coord, "result " + newMdl, clid);

    }

    public static double sigmoid(double x, double mu, double stddev, double a, double b) {
        return (a - b) * (1 - normp((x - mu) / stddev)) + b;
    }

    public static double log2(double x) {
        final double l2 = 1.0 / Math.log(2.0);
        return Math.log(x) * l2;
    }

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }

    public static double normp(double z) {

        double zabs;
        double p;
        double expntl, pdf;

        final double p0 = 220.2068679123761;
        final double p1 = 221.2135961699311;
        final double p2 = 112.0792914978709;
        final double p3 = 33.91286607838300;
        final double p4 = 6.373962203531650;
        final double p5 = .7003830644436881;
        final double p6 = .3526249659989109E-01;

        final double q0 = 440.4137358247522;
        final double q1 = 793.8265125199484;
        final double q2 = 637.3336333788311;
        final double q3 = 296.5642487796737;
        final double q4 = 86.78073220294608;
        final double q5 = 16.06417757920695;
        final double q6 = 1.755667163182642;
        final double q7 = .8838834764831844E-1;

        final double cutoff = 7.071;
        final double root2pi = 2.506628274631001;

        zabs = Math.abs(z);
        if (z > 37.0) {
            p = 1.0;
            return p;
        }
        if (z < -37.0) {
            p = 0.0;
            return p;
        }
        expntl = Math.exp(-.5 * zabs * zabs);
        pdf = expntl / root2pi;
        if (zabs < cutoff) {
            p = expntl * ((((((p6 * zabs + p5) * zabs + p4) * zabs + p3) * zabs
                    + p2) * zabs + p1) * zabs + p0) / (((((((q7 * zabs + q6) * zabs
                    + q5) * zabs + q4) * zabs + q3) * zabs + q2) * zabs + q1) * zabs
                    + q0);
        } else {
            p = pdf / (zabs + 1.0 / (zabs + 2.0 / (zabs + 3.0 / (zabs + 4.0
                    / (zabs + 0.65)))));
        }
        if (z < 0.0) {
            return p;
        } else {
            p = 1.0 - p;
            return p;
        }
    }
}
