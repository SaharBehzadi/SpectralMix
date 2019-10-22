/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import java.util.Random;
import Jama.*;

public class NewMainChris {

    public static class pair implements Comparable {

        public double dist;
        public boolean isEdge;

        public pair(double dist, boolean isEdge) {
            this.dist = dist;
            this.isEdge = isEdge;
        }

        public int compareTo(Object o) {
            pair other = (pair) o;
            if (this.dist < other.dist) {
                return -1;
            }
            if (this.dist > other.dist) {
                return 1;
            }
            return 0;
        }
    }

    public static void main(String[] args) {
        test_sigmoid();
        Random r = new Random(1);
        int n = 300;
        int d = 2;
        int k = 4;
        int histbins = 30;
        double[][] db = new double[n][d];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < d; j++) {
                db[i][j] = r.nextDouble();
            }
        }
        boolean[][] adjacency = new boolean[n][n];
        for (int i = 0; i < n; i++) {
            for (int kk = 0; kk < k; kk++) {
                double mindist = Double.MAX_VALUE;
                int minj = -1;
                for (int j = 0; j < n; j++) {
                    if (i != j && !adjacency[i][j]) {
                        double actdist = dist(db[i], db[j]);
                        if (actdist < mindist) {
                            minj = j;
                            mindist = actdist;
                        }

                    }
                }
                adjacency[i][minj] = true;
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                adjacency[i][j] |= adjacency[j][i];
            }
        }
        for (int i = 0; i < n; i++) {
            System.out.print(adjacency[i][0] ? 1 : 0);
            for (int j = 1; j < n; j++) {
                System.out.print("," + (adjacency[i][j] ? 1 : 0));
            }
            System.out.println();
        }
        int count = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (adjacency[i][j]) {
                    count++;
                }
            }
        }
        System.out.println(count + "/" + (n * (n - 1) / 2) + " Links");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (adjacency[i][j]) {
                    System.out.println(db[i][0] + "," + db[i][1]);
                    System.out.println(db[j][0] + "," + db[j][1]);
                    System.out.println();
                }
            }
        }
        for (int i = 0; i < n; i++) {
            System.out.println(db[i][0] + "," + db[i][1]);
        }
        double[][] dist = floydWarshall(adjacency);
        for (int i = 0; i < n; i++) {
            System.out.print((int) dist[i][0]);
            for (int j = 1; j < n; j++) {
                System.out.print("," + (int) dist[i][j]);
            }
            System.out.println();
        }

        // TEST TEST TEST Use original distances
//        for (int i=0 ; i<n ; i++)
//            for (int j=0 ; j<n ; j++)
//                dist[i][j]=dist(db[i], db[j]) ;
        // END TEST

        // double centering
        double[] colsum = new double[n];
        double allsum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                colsum[i] += dist[i][j] * dist[i][j];
            }
            allsum += colsum[i];
        }
        double[][] bmatrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                bmatrix[i][j] = -0.5 * (dist[i][j] * dist[i][j] - (colsum[i] + colsum[j]) / n + allsum / n / n);
            }
        }

        EigenvalueDecomposition evd = new EigenvalueDecomposition(new Matrix(bmatrix));
        int bestEV = 0;
        for (int i = 1; i < n; i++) {
            if (evd.getD().get(i, i) > evd.getD().get(bestEV, bestEV)) {
                bestEV = i;
            }
        }
        int secondEV = bestEV == 0 ? 1 : 0;
        for (int i = secondEV + 1; i < n; i++) {
            if (i != bestEV && evd.getD().get(i, i) > evd.getD().get(secondEV, secondEV)) {
                secondEV = i;
            }
        }
        Matrix v = evd.getV();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (adjacency[i][j]) {
//                    System.out.println(v.get(n-1, i)+","+v.get(n-2, i)); PROVEN WRONG
//                    System.out.println(v.get(n-1, j)+","+v.get(n-2, j)); PROVEN WRONG
                    System.out.println(v.get(i, bestEV) + "," + v.get(i, secondEV));
                    System.out.println(v.get(j, bestEV) + "," + v.get(j, secondEV));
                    System.out.println();
                }
            }
        }
        double[][] dbRestored = new double[n][d];
        for (int i = 0; i < n; i++) {
            dbRestored[i][0] = v.get(i, bestEV) * Math.sqrt(evd.getD().get(bestEV, bestEV));
            dbRestored[i][1] = v.get(i, secondEV) * Math.sqrt(evd.getD().get(secondEV, secondEV));
        }



        // TEST TEST TEST
        // DETERMINE MAJORIZATION


        majorization(dbRestored, adjacency, 1000);
        //weightedMajorization2 (dbRestored, adjacency, 10000) ;


        System.out.println("After Majorization");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (adjacency[i][j]) {
                    System.out.println(dbRestored[i][0] + "," + dbRestored[i][1]);
                    System.out.println(dbRestored[j][0] + "," + dbRestored[j][1]);
                    System.out.println();
                }
            }
        }

        // END TEST MAJORIZATION

         //TEST TEST TEST
         //random init
//        Random rr = new Random (216) ;
//        for (int i=0 ; i<n ; i++)
//            for (int j=0 ; j<d ; j++)
//                dbRestored[i][j] = rr.nextGaussian() ;
         // END TEST

        improveCoordinatesMedian(dbRestored, adjacency);
        System.out.println("After Improve");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (adjacency[i][j]) {
                    System.out.println(dbRestored[i][0] + "," + dbRestored[i][1]);
                    System.out.println(dbRestored[j][0] + "," + dbRestored[j][1]);
                    System.out.println();
                }
            }
        }

        double minx, miny, maxx, maxy;
        minx = maxx = dbRestored[0][0];
        miny = maxy = dbRestored[0][1];
        for (int i = 1; i < n; i++) {
            minx = Math.min(minx, dbRestored[i][0]);
            maxx = Math.max(maxx, dbRestored[i][0]);
            miny = Math.min(miny, dbRestored[i][1]);
            maxy = Math.max(maxy, dbRestored[i][1]);
        }
        maxx *= 1.0 + Math.pow(2.0, -50);
        maxy *= 1.0 + Math.pow(2.0, -50);
        double mindist = dist(dbRestored[0], dbRestored[1]);
        double maxdist = mindist;
        for (int i = 2; i < n; i++) {
            for (int j = 0; j < i; j++) {
                mindist = Math.min(mindist, dist(dbRestored[i], dbRestored[j]));
                maxdist = Math.max(maxdist, dist(dbRestored[i], dbRestored[j]));
            }
        }
        maxdist *= 1.0 + Math.pow(2.0, -50);
        int[] histEdge = new int[histbins];
        int[] histAll = new int[histbins];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                int bin = (int) ((dist(dbRestored[i], dbRestored[j]) - mindist) / (maxdist - mindist) * histbins);
                histAll[bin]++;
                if (adjacency[i][j]) {
                    histEdge[bin]++;
                }
            }
        }
        int histAllAll = histAll[0];
        int histEdgeAll = histEdge[0];
        for (int i = 1; i < histbins; i++) {
            histAllAll += histAll[i];
            histEdgeAll += histEdge[i];
        }
        if (histAllAll != n * (n - 1) / 2) {
            System.out.println("ERROR ERROR ERROR " + histAllAll);
        }
        double costAll = 0;
        for (int i = 0; i < histbins; i++) {
            costAll += entr(histEdge[i], histAll[i]) * histAll[i] / histAllAll;
        }

        System.out.println("HISTOGRAM:\n-------------");
        for (int i = 0; i < histbins; i++) {
            System.out.println(histEdge[i] + ", " + histAll[i]);
        }

        // Order the distances
        pair[] pairs = new pair[n * (n - 1) / 2];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                pairs[i * (i - 1) / 2 + j] = new pair(dist(dbRestored[i], dbRestored[j]), adjacency[i][j]);
            }
        }
        java.util.Arrays.sort(pairs);
        double sum = 0;
        double sqrsum = 0;
        double weight = 0;
        for (int i = 1; i < n * (n - 1) / 2; i++) {
            if (pairs[i - 1].isEdge && !pairs[i].isEdge) {
                weight += 1;
                double avgdist = (pairs[i - 1].dist + pairs[i].dist) / 2;
                sum += avgdist;
                sqrsum += avgdist * avgdist;
            } else if (!pairs[i - 1].isEdge && pairs[i].isEdge) {
                weight -= 1;
                double avgdist = (pairs[i - 1].dist + pairs[i].dist) / 2;
                sum -= avgdist;
                sqrsum -= avgdist * avgdist;
            }
        }
        double stddev = Math.sqrt(sqrsum / weight - sum * sum / weight * weight);
        System.out.println(sum / weight + " +- " + stddev + " (sum: " + sum + "; weight: " + weight + "; sqrsum: "
                + sqrsum + "; mindist: " + mindist + "; maxdist: " + maxdist + ")");

        // and now determine the scaling parameters a and b of the sigmoid

        double sigavg = 0.0;
        double edgeprob = (double) histEdgeAll / histAllAll;
        for (int i = 0; i < n * (n - 1) / 2; i++) {
            sigavg += sigmoid(pairs[i].dist, sum / weight, stddev, 1, 0);
        }
        sigavg /= n * (n - 1) / 2;
        double nominator = 0;
        double denominator = 0;
        for (int i = 0; i < n * (n - 1) / 2; i++) {
            double h = sigmoid(pairs[i].dist, sum / weight, stddev, 1, 0) - sigavg;
            nominator += h * ((pairs[i].isEdge ? 1 : 0) - edgeprob);
            denominator += h * h;
        }
        double b = edgeprob - nominator / denominator * sigavg;
        double a = nominator / denominator + b;
        System.out.println("a=" + a + "; b=" + b);

        System.out.println("\nHISTOGRAM:\n-------------");
        for (int i = 0; i < histbins; i++) {
            double h = ((double) i + 0.5) / histbins * (maxdist - mindist) + mindist;
            System.out.println(histEdge[i] + ", " + histAll[i] + ", "
                    + (double) histEdge[i] / histAll[i] + ", " + sigmoid(h, sum / weight, stddev, a, b));
        }

                    PairOwn[] p = new PairOwn[n*(n-1)/2];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(dbRestored[i], dbRestored[j]), adjacency[i][j]);
                }
            }
            Sigmoid s = new Sigmoid(p);
            double sigcost = 0 ;
            for (int i=0 ; i<n*(n-1)/2 ; i++)
                    if (p[i].isEdge)
                        sigcost -= log2(s.f(p[i].dist)) ;
                    else
                        sigcost -= log2(1-s.f(p[i].dist)) ;


        System.out.println("Overall Entropy of reconstructed DB: " + costAll);
        System.out.println(" - Bitvolume of Adj Matrix: " + costAll * histAllAll);
        System.out.println(" - Coding of coord acc. to BIC: " + n * Math.log(histAllAll / n) / Math.log(2));
        System.out.println(" - Bitvolume together: " + (costAll * histAllAll + n * Math.log(histAllAll / n) / Math.log(2)));

        // Cross-Check Coding Cost Using the Sigmoid Function

        double sigmoidCost = 0.0;
        for (int i = 0; i < n * (n - 1) / 2; i++) {
            if (pairs[i].isEdge) {
                sigmoidCost -= Math.log(sigmoid(pairs[i].dist, sum / weight, stddev, a, b)) / Math.log(2);
            } else {
                sigmoidCost -= Math.log(1 - sigmoid(pairs[i].dist, sum / weight, stddev, a, b)) / Math.log(2);
            }
        }
        System.out.println(" - Bitvolume of Adj. Matrix with SIGMOID: " + sigmoidCost);
        System.out.println("\nEntropy of " + histEdgeAll + " edges in " + histAllAll + " possible Edges: " + entr(histEdgeAll, histAllAll));
        System.out.println(" - Bitvolume of Adj Matrix: " + entr(histEdgeAll, histAllAll) * histAllAll);



        // Entropy in original coordinates:
        mindist = dist(db[0], db[1]);
        maxdist = mindist;
        for (int i = 2; i < n; i++) {
            for (int j = 0; j < i; j++) {
                mindist = Math.min(mindist, dist(db[i], db[j]));
                maxdist = Math.max(maxdist, dist(db[i], db[j]));
            }
        }
        maxdist *= 1.0 + Math.pow(2.0, -50);
        histEdge = new int[histbins];
        histAll = new int[histbins];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                int bin = (int) ((dist(db[i], db[j]) - mindist) / (maxdist - mindist) * histbins);
                histAll[bin]++;
                if (adjacency[i][j]) {
                    histEdge[bin]++;
                }
            }
        }
        histAllAll = histAll[0];
        histEdgeAll = histEdge[0];
        for (int i = 1; i < histbins; i++) {
            histAllAll += histAll[i];
            histEdgeAll += histEdge[i];
        }
        if (histAllAll != n * (n - 1) / 2) {
            System.out.println("ERROR ERROR ERROR " + histAllAll);
        }
        costAll = 0;
        for (int i = 0; i < histbins; i++) {
            costAll += entr(histEdge[i], histAll[i]) * histAll[i] / histAllAll;
        }
        System.out.println("Overall Entropy of original coordinates: " + costAll);
        System.out.println(" - Bitvolume of Adj Matrix: " + costAll * histAllAll);

        // Entropy in quantized restored coordinates

        for (int bits = 10; bits > 0; bits--) {
            // truncate coordinates to ... bits
            for (int i = 0; i < n; i++) {
                dbRestored[i][0] = ((int) ((dbRestored[i][0] - minx) / (maxx - minx) * Math.pow(2.0, bits)))
                        / Math.pow(2.0, bits) * (maxx - minx) + minx;
                dbRestored[i][1] = ((int) ((dbRestored[i][1] - miny) / (maxy - miny) * Math.pow(2.0, bits)))
                        / Math.pow(2.0, bits) * (maxy - miny) + miny;
            }
            mindist = dist(dbRestored[0], dbRestored[1]);
            maxdist = mindist;
            for (int i = 2; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    mindist = Math.min(mindist, dist(dbRestored[i], dbRestored[j]));
                    maxdist = Math.max(maxdist, dist(dbRestored[i], dbRestored[j]));
                }
            }
            maxdist *= 1.0 + Math.pow(2.0, -50);
            histEdge = new int[histbins];
            histAll = new int[histbins];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    int bin = (int) ((dist(dbRestored[i], dbRestored[j]) - mindist) / (maxdist - mindist) * histbins);
                    histAll[bin]++;
                    if (adjacency[i][j]) {
                        histEdge[bin]++;
                    }
                }
            }
            histAllAll = histAll[0];
            histEdgeAll = histEdge[0];
            for (int i = 1; i < histbins; i++) {
                histAllAll += histAll[i];
                histEdgeAll += histEdge[i];
            }
            if (histAllAll != n * (n - 1) / 2) {
                System.out.println("ERROR ERROR ERROR " + histAllAll);
            }
            costAll = 0;
            for (int i = 0; i < histbins; i++) {
                costAll += entr(histEdge[i], histAll[i]) * histAll[i] / histAllAll;
            }
            System.out.println("Overall Entropy of " + bits + " bit quantized coordinates: " + costAll);
            System.out.println(" - Bitvolume of Adj Matrix: " + costAll * histAllAll);
            System.out.println(" - Bitvolume of coordinates: " + 2 * n * bits);
            System.out.println(" - Bitvolume of both: " + (2 * n * bits + costAll * histAllAll));
        }

        // Entropy in quantized original coordinates
        minx = maxx = db[0][0];
        miny = maxy = db[0][1];
        for (int i = 1; i < n; i++) {
            minx = Math.min(minx, db[i][0]);
            maxx = Math.max(maxx, db[i][0]);
            miny = Math.min(miny, db[i][1]);
            maxy = Math.max(maxy, db[i][1]);
        }
        maxx *= 1.0 + Math.pow(2.0, -50);
        maxy *= 1.0 + Math.pow(2.0, -50);

        for (int bits = 10; bits > 0; bits--) {
            // truncate coordinates to ... bits
            for (int i = 0; i < n; i++) {
                db[i][0] = ((int) ((db[i][0] - minx) / (maxx - minx) * Math.pow(2.0, bits)))
                        / Math.pow(2.0, bits) * (maxx - minx) + minx;
                db[i][1] = ((int) ((db[i][1] - miny) / (maxy - miny) * Math.pow(2.0, bits)))
                        / Math.pow(2.0, bits) * (maxy - miny) + miny;
            }
            mindist = dist(db[0], db[1]);
            maxdist = mindist;
            for (int i = 2; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    mindist = Math.min(mindist, dist(db[i], db[j]));
                    maxdist = Math.max(maxdist, dist(db[i], db[j]));
                }
            }
            maxdist *= 1.0 + Math.pow(2.0, -50);
            histEdge = new int[histbins];
            histAll = new int[histbins];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    int bin = (int) ((dist(db[i], db[j]) - mindist) / (maxdist - mindist) * histbins);
                    histAll[bin]++;
                    if (adjacency[i][j]) {
                        histEdge[bin]++;
                    }
                }
            }
            histAllAll = histAll[0];
            histEdgeAll = histEdge[0];
            for (int i = 1; i < histbins; i++) {
                histAllAll += histAll[i];
                histEdgeAll += histEdge[i];
            }
            if (histAllAll != n * (n - 1) / 2) {
                System.out.println("ERROR ERROR ERROR " + histAllAll);
            }
            costAll = 0;
            for (int i = 0; i < histbins; i++) {
                costAll += entr(histEdge[i], histAll[i]) * histAll[i] / histAllAll;
            }
            System.out.println("Overall Entropy of " + bits + " bit original coordinates: " + costAll);
            System.out.println(" - Bitvolume of Adj Matrix: " + costAll * histAllAll);
            System.out.println(" - Bitvolume of coordinates: " + 2 * n * bits);
            System.out.println(" - Bitvolume of both: " + (2 * n * bits + costAll * histAllAll));
        }
            System.out.println("SIGMOID: "+s) ;
            System.out.println ("SIGCOST: "+sigcost) ;

    }

    public static void test_sigmoid() {
        int n = 1000000;
        int histbins = 100;
        // Parameters of generating Sigmoid function
        Sigmoid gen = new Sigmoid(1.0, 1.0, 0.8, 0.2);
        double minG = 0.0;
        double maxG = 10.0;
        PairOwn[] pairs = new PairOwn[n];
        for (int i = 0; i < n; i++) {
            pairs[i] = gen.generateRandom(minG, maxG);
        }
        Sigmoid detected = new Sigmoid(pairs);
        System.out.println(gen + " =?= " + detected);
        int[] histEdge = new int[histbins];
        int[] histAll = new int[histbins];
        for (int i = 1; i < n; i++) {
            int bin = (int) ((pairs[i].dist - minG) / (maxG - minG) * histbins);
            histAll[bin]++;
            if (pairs[i].isEdge) {
                histEdge[bin]++;
            }
        }
        for (int i = 0; i < histbins - 1; i++) {
            double x = (i + 0.5) * (maxG - minG) / histbins;
            System.out.println((double) histEdge[i] / histAll[i] + ", " + gen.f(x) + ", " + (-(double) histEdge[i + 1] / histAll[i + 1] + (double) histEdge[i] / histAll[i]));
        }
        double sum = 0.0;
        double sumsq = 0.0;
        double weight = 0.0;
        for (int i = 0; i < histbins - 1; i++) {
            double x = (i + 1.0) * (maxG - minG) / histbins;
            double h = (double) histEdge[i] / histAll[i] - (double) histEdge[i + 1] / histAll[i + 1];
            weight += h;
            sum += h * x;
            sumsq += h * x * x;
        }
        System.out.println("w=" + weight + "; mu=" + sum / weight + "; sigma=" + (sumsq / weight - sum * sum / weight / weight));
    }

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }

    public static double entr(int pos, int all) {
        if (pos == 0 || pos == all) {
            return 0.0;
        }
        return -(Math.log((double) pos / all) * ((double) pos / all)
                + Math.log((double) (all - pos) / all) * (double) (all - pos) / all)
                / Math.log(2.0);
    }

    public static double[][] floydWarshall(boolean[][] adjacency) {
        int d = adjacency.length;
        double[][] result = new double[d][d];
        for (int i = 1; i < d; i++) {
            for (int j = 0; j < i; j++) {
                result[i][j] = result[j][i] = adjacency[i][j] ? 1 : 30;
            }
        }
        for (int k = 0; k < d; k++) {
            for (int i = 0; i < d; i++) {
                for (int j = 0; j < d; j++) {
                    result[i][j] = Math.min(result[i][j], result[i][k] + result[k][j]);
                }
            }
        }

        return result;
    }

    public static void improveCoordinates(double[][] coord, boolean[][] adjacency) {
        int n = coord.length;
        int d = coord[0].length;
        int m = n * (n - 1) / 2;
        Random r = new Random(2);
        if (adjacency.length != n || adjacency[0].length != n) {
            throw (new IndexOutOfBoundsException());
        }
        for (double sig = 0.05; sig > 0.000001; sig *= 0.995) {
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), adjacency[i][j]);
                }
            }
            Sigmoid s = new Sigmoid(p);
            for (int i = 0; i < n; i++) {
                // now determine individual cost of all edges and non-edges of node i
                double[] newCoord = coord[i].clone();
                double costBefore = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        if (adjacency[i][j]) {
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
                        if (adjacency[i][j]) {
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
    }

    public static void improveCoordinatesSteepestDescent(double[][] coord, boolean[][] adjacency) {
        int n = coord.length;
        int d = coord[0].length;
        int m = n * (n - 1) / 2;
        Random r = new Random(2);
        if (adjacency.length != n || adjacency[0].length != n) {
            throw (new IndexOutOfBoundsException());
        }
        for (int loop = 0 ; loop<1000 ; loop++) {
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), adjacency[i][j]);
                }
            }
            Sigmoid s = new Sigmoid(p);
            for (int i = 0; i < n; i++) {
                // now determine individual cost of all edges and non-edges of node i
                double[] newCoord = new double[d];
                double sumWeight = 0.0 ;
                double costBefore = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        if (adjacency[i][j]) {
                            costBefore -= log2(s.f(dist(coord[i], coord[j])));
                        } else {
                            costBefore -= log2(1.0 - s.f(dist(coord[j], coord[j])));
                        }
                    }
                }
                for (int j = 0; j < n; j++) if(i!=j){
                    double weight = 0;
                    if (adjacency[i][j])
                        weight = s.diffCostEdge(dist(coord[i],coord[j])) ;
                    else
                        weight = s.diffCostNoEdge(dist(coord[i],coord[j])) ;
                    // weight/= dist(coord[i],coord[j]) ;
                    if (weight<0){
                        weight = -weight ;
                        for (int jj=0 ; jj<d ; jj++)
                            newCoord[jj] += weight * (2*coord[i][jj]-coord[j][jj]) ;
                    } else
                    for (int jj=0; jj<d ; jj++)
                        newCoord[jj] += weight * coord[j][jj] ;
                    sumWeight += weight ;
                }
                // now determine new cost
                double costAfter = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        if (adjacency[i][j]) {
                            costAfter -= log2(s.f(dist(newCoord, coord[j])));
                        } else {
                            costAfter -= log2(1.0 - s.f(dist(newCoord, coord[j])));
                        }
                    }
                }
                if (costAfter > costBefore) {
                    loop=loop; // place for breakpoint
                }
                for (int jj=0 ; jj<d ; jj++)
                    coord[i][jj] = 0.95*coord[i][jj] +0.05*newCoord[jj]/sumWeight ;
                //}
            }
        }
    }

    public static void improveCoordinatesMedian(double[][] coord, boolean[][] adjacency) {
        int n = coord.length;
        int d = coord[0].length;
        int m = n * (n - 1) / 2;
        if (adjacency.length != n || adjacency[0].length != n) {
            throw (new IndexOutOfBoundsException());
        }
        for (int loop = 0 ; loop<1000 ; loop++) {
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), adjacency[i][j]);
                }
            }
            Sigmoid s = new Sigmoid(p);
            for (int i = 0; i < n; i++) {
                // now determine individual cost of all edges and non-edges of node i
                double sumWeight = 0.0 ;
                PairForMedian [][] pfm = new PairForMedian[d][n-1] ;
                for (int j = 0; j < n; j++) if(i!=j){
                    double weight = 0;
                    if (adjacency[i][j])
                        weight = s.diffCostEdge(dist(coord[i],coord[j])) ;
                    else
                        weight = s.diffCostNoEdge(dist(coord[i],coord[j])) ;
                    // weight/= dist(coord[i],coord[j]) ;
                    if (weight<0){
                        weight = -weight ;
                        for (int jj=0 ; jj<d ; jj++)
                            pfm[jj][j-(j>i?1:0)] = new PairForMedian(2*coord[i][jj]-coord[j][jj],weight) ;
                    } else
                    for (int jj=0; jj<d ; jj++)
                        pfm[jj][j-(j>i?1:0)] = new PairForMedian(coord[j][jj],weight) ;
                    sumWeight += weight ;
                }
                // now determine new cost
                for (int jj=0 ; jj<d ; jj++) {
                     java.util.Arrays.sort(pfm[jj]) ;
                     double sumHalfWeight = 0;
                     int j=0 ;
                     for (j=0 ; j<n-2&&sumHalfWeight<sumWeight/2 ; j++)
                         sumHalfWeight += pfm[jj][j].weight ;
                     coord[i][jj] = 0.9 * coord[i][jj] +0.1 * pfm[jj][j].value ;
                }
            }
        }
    }

    public static void majorization(double[][] dbRestored, boolean[][] adjacency, int iterations) {
        int n = dbRestored.length;
        int d = dbRestored[0].length;
        //int m = n*(n-1)/2 ;
        for (int loop = 0; loop < iterations; loop++) {

            double[][] dbRestoredNew = new double[n][d];
            double[][] sij = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    sij[i][j] = 0;
                    double h = dist(dbRestored[i], dbRestored[j]);
                    if (h != 0) {
                        sij[i][j] = (adjacency[i][j] ? 0 : 1) / h;
                    }
                }
            }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        for (int jj = 0; jj < d; jj++) {
                            dbRestoredNew[i][jj] += dbRestored[j][jj] + sij[i][j] * (dbRestored[i][jj] - dbRestored[j][jj]);
                        }
                    }
                }
            }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    dbRestored[i][j] = dbRestoredNew[i][j] / (n - 1);
                }
            }

        }
    }

    public static void weightedMajorization(double[][] dbRestored, boolean[][] adjacency, int iterations) {
        int n = dbRestored.length;
        int d = dbRestored[0].length;
        int m = n * (n - 1) / 2;
        int bestiter = -1 ;
        double bestcost = Double.MAX_VALUE ;
        for (int loop = 0; loop < iterations; loop++) {
            double sumSumWeight = 0;
            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(dbRestored[i], dbRestored[j]), adjacency[i][j]);
                }
            }
            Sigmoid s = new Sigmoid(p);

            double[][] dbRestoredNew = new double[n][d];
            double[][] sij = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    sij[i][j] = 0;
                    double h = dist(dbRestored[i], dbRestored[j]);
                    if (h != 0) {
                        sij[i][j] = (adjacency[i][j] ? 0 : 1) / h;
                    }
                }
            }
            for (int i = 0; i < n; i++) {
                double sumWeight = 0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        double h=0 ;
                        if (adjacency[i][j]) {
                            h = -log2(s.f(dist(dbRestored[i], dbRestored[j])));
                        } else {
                            h = -log2(1.0-s.f(dist(dbRestored[i], dbRestored[j])));
                        }
                        sumWeight += h;
                        for (int jj = 0; jj < d; jj++) {
                            dbRestoredNew[i][jj] += h * (dbRestored[j][jj] + sij[i][j] * (dbRestored[i][jj] - dbRestored[j][jj]));
                        }
                    }
                }
                for (int jj = 0; jj < d; jj++) {
                    dbRestoredNew[i][jj] /= sumWeight;
                }
                sumSumWeight += sumWeight ;
            }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    dbRestored[i][j] = 0.9*dbRestored[i][j]+0.1*dbRestoredNew[i][j] ;
                }
            }
            if(loop % (iterations/10)==0)
                System.out.println ("Iteration No "+loop+" / "+iterations+": "+sumSumWeight+ "(Best: "+bestcost+" after "+bestiter+")");
            if (sumSumWeight < bestcost) {
                bestcost = sumSumWeight ;
                bestiter = loop ;
            }

        }
    }

    public static void weightedMajorization2(double[][] dbRestored, boolean[][] adjacency, int iterations) {
        int n = dbRestored.length;
        int d = dbRestored[0].length;
        for (int loop = 0; loop < iterations; loop++) {
            double[][] dbRestoredNew = new double[n][d];
            double[][] sij = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    sij[i][j] = 0;
                    double h = dist(dbRestored[i], dbRestored[j]);
                    if (h != 0) {
                        sij[i][j] = 1 / h;
                    }
                }
            }
            for (int i = 0; i < n; i++) {
                double sumWeight = 0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        double h=0 ;
                        if (adjacency[i][j]) {
                            h = 1 ;
                        } else {
                            h = -1 ;
                        }
                        sumWeight += h;
                        for (int jj = 0; jj < d; jj++) {
                            dbRestoredNew[i][jj] += h * (dbRestored[j][jj] + sij[i][j] * (dbRestored[i][jj] - dbRestored[j][jj]));
                        }
                    }
                }
                for (int jj = 0; jj < d; jj++) {
                    dbRestoredNew[i][jj] /= sumWeight;
                }
            }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < d; j++) {
                    dbRestored[i][j] = dbRestoredNew[i][j] ;
                }
            }
            if(loop % (iterations/10)==0)
                System.out.println ("Iteration No "+loop+" / "+iterations);
        }
    }

    public static double sigmoid(double x, double mu, double stddev, double a, double b) {
        return (a - b) * (1 - normp((x - mu) / stddev)) + b;
    }

    public static double log2(double x) {
        final double l2 = 1.0 / Math.log(2.0);
        return Math.log(x) * l2;
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

