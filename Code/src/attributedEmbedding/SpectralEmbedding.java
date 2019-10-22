/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package attributedEmbedding;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import edu.uci.ics.jung.graph.Graph;
import java.util.Collection;
import java.util.Random;
import nature.IO;
import nature.Visualization;

/**
 *
 * @author plantc59cs
 */
public class SpectralEmbedding {

    Graph g;
    int d; //embedding dimensionality
    double[][] objCoord; //numVertices x d
    int[] classId;
    int[] sumWeights; //for each object: the number of neighbors plus the number of categories; missing values possible
    int numObj;
    static boolean writeResult = true;

//    public JointEmbedding(Graph g, int[][] attributes, int[] numCat, int d) {
//        this.g = g;
//        this.attributes = attributes;
//        this.numCat = numCat;
//        this.d = d;
//    }
    public SpectralEmbedding(Graph g, int d, int[] classId) {
        this.g = g;
        this.d = d;
        this.classId = classId;
        numObj = g.getVertexCount();
        objCoord = new double[numObj][d];
        sumWeights = new int[numObj];
        for (int i = 0; i < numObj; i++) {
            sumWeights[i] += g.getNeighborCount(i);

        }

    }

    public void checkEigenvector() {
        //random walk normalization
        double[][] lap = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = i + 1; j < numObj; j++) {
                if (g.isNeighbor(i, j)) {
                    lap[i][j] = lap[j][i] = -1;
                }
            }
            lap[i][i] = g.degree(i);
        }
        //laplacian = new Matrix(lap);
        double[][] dm = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            dm[i][i] = 1.0 / g.degree(i);
           // dm[i][i] = lap[i][i];
        }
        Matrix laplacian = new Matrix(dm).times(new Matrix(lap));
        EigenvalueDecomposition eig = new EigenvalueDecomposition(laplacian);
        double[][] ew = eig.getD().getArrayCopy();
        double[][] ev = eig.getV().getArrayCopy();
        if (writeResult) {
            IO ea = new IO();
            ea.writeDoubleToMatlab(ew, "ew", "ew");
            ea.writeDoubleToMatlab(ev, "ev", "ev");

        }
    }

    //init object coordinates randomly and category coords at the center of all assigned objects
    public void init(int seed) {
        Random r = new Random(seed);
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < d; j++) {
                objCoord[i][j] = r.nextDouble();
            }
        }

        displayGraphClass();
    }

    private void displayGraphClass() {
        Visualization v = new Visualization(g);
        double[][] cc = new double[numObj][2];
        for (int i = 0; i < numObj; i++) {
            cc[i][0] = objCoord[i][1];
            cc[i][1] = objCoord[i][2];
        }
        v.displayCoordNew(cc, "Spectral Class", classId);
    }

    private double obj() {
        double graphCost = 0.0;
        for (int i = 0; i < numObj; i++) {
            for (int j = i + 1; j < numObj; j++) {
                if (g.isNeighbor(i, j)) {
                    double dist = 0.0;
                    for (int l = 0; l < d; l++) {
                        dist += (objCoord[i][l] - objCoord[j][l]) * (objCoord[i][l] - objCoord[j][l]);
                    }
                    //unweightedGraph
                    graphCost += dist;
                }
            }
        }

        double cost = graphCost;
        System.out.println(cost);
        return cost;
    }

    public void run() {
        boolean converged = false;
        double minCost = obj();
        int iter = 0;
        while (!converged || iter < 5000) {
            updateObjects();

            double aktCost = obj();
            if (aktCost < minCost) {
                minCost = aktCost;
            } else {
                converged = true;
            }
            iter++;
        }
        if (writeResult) {
            IO ea = new IO();
            ea.writeDoubleToMatlab(objCoord, "objCoord", "objCoord");

        }

        displayGraphClass();

    }

    private void updateObjects() {
        double[][] coordNew = new double[numObj][d];
        for (int i = 0; i < numObj; i++) {
            //graph contribution
            Collection<Integer> n = g.getNeighbors(i);
            for (int p : n) {
                for (int j = 0; j < d; j++) {
                    if (p != i) {
                        coordNew[i][j] += objCoord[p][j] / sumWeights[i];
                    }
                }
            }
            //attribute contribution
//            for(int k = 0; k < numAtt; k++){
//                for(int j = 0; j < d; j++){
//                    coordNew[i][j] += catCoord[startIndex[k] + attributes[i][k]][j]/sumWeights[i];
//                }
//            }

        }
        //TEST
        //double bla = obj(coordNew);
        //TEST
        objCoord = modifiedGramSchmidt(coordNew);
        //coord = coordNew;
        //TEST
        //double blabla = obj(coord);
        //TEST

    }

    //orthogonalize the column vectors in v
    private double[][] modifiedGramSchmidt(double[][] v) {
        int k = v[0].length;
        int vLength = v.length;
        for (int j = 0; j < k; j++) {
            for (int i = 0; i < j; i++) {
                double skalarprod = 0.0;
                double self_i = 0.0;
                double proj_vi_vj = 0.0;
                for (int l = 0; l < vLength; l++) {
                    skalarprod += v[l][i] * v[l][j];
                    self_i += v[l][i] * v[l][i];
                }
                for (int l = 0; l < vLength; l++) {
                    proj_vi_vj = (skalarprod / self_i) * v[l][i];
                    v[l][j] = v[l][j] - proj_vi_vj;
                }
            } //i
            double norm_j = 0.0;
            for (int l = 0; l < vLength; l++) {
                norm_j += v[l][j] * v[l][j];
            }
            norm_j = Math.sqrt(norm_j);
            for (int l = 0; l < vLength; l++) {
                v[l][j] = v[l][j] / norm_j;
            }
        }//j

        // writeTest(v, "orth.mat");
        return v;
    }

}
