/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package attributedEmbedding;

import datautils.MyEdge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import java.util.Collection;
import java.util.Random;
import nature.IO;
import nature.Visualization;

/**
 *
 * @author plantc59cs
 */
public class JointEmbedding {

    Graph g;
    int[][] attributes; //n x numAtt; each categorical attribute value encoded as integer starting from 0
//    int[] singleAtt; //n; for the case of only one categorical attribute
//    int singleCat; //number of categories of the only categorical attribute
    int numAtt; //number of categorical attributes 
    //int[] numCat; //number of Categories of each attribute
    CountPack[] catCount;
    double[][] catCoord; //numCat x d: coords of category centers
    int d; //embedding dimensionality
    double[][] objCoord; //numVertices x d
    int[] classId;
    int[] startIndex; //starting index of each Attribute in catCoord
    double[] sumWeights; //for each object: the number of neighbors plus the number of categories; missing values possible
    int numObj;
    double weightGraph;
    double weightAttributes;
    static boolean writeResult = true;
    static boolean weighted = true;

//    public JointEmbedding(Graph g, int[][] attributes, int[] numCat, int d) {
//        this.g = g;
//        this.attributes = attributes;
//        this.numCat = numCat;
//        this.d = d;
//    }
    public JointEmbedding(Graph g, int[] singleAtt, int singleCat, int d, int[] classId) {
        this.g = g;
        this.d = d;
        this.classId = classId;
        numObj = g.getVertexCount();
        numAtt = 1;
        catCount = new CountPack[1];
        attributes = new int[numObj][1];
        int[] counts = new int[singleCat];
        for (int i = 0; i < numObj; i++) {
            attributes[i][0] = singleAtt[i];
            counts[singleAtt[i]]++;
        }
        catCount[0] = new CountPack(counts);
        objCoord = new double[numObj][d];
        catCoord = new double[singleCat][d];
        startIndex = new int[1];
        startIndex[0] = 0;
        sumWeights = new double[numObj];
        for (int i = 0; i < numObj; i++) {
            sumWeights[i] += g.getNeighborCount(i);
            sumWeights[i]++; //for the categorical attribute
        }

    }

    public JointEmbedding(Graph g, int[] singleAtt, int singleCat, int d, int[] classId, double weightGraph, double weightAttributes) {
        this.g = g;
        this.d = d;
        this.classId = classId;
        this.weightGraph = weightGraph;
        this.weightAttributes = weightAttributes;
        System.out.println("weightGraph: " + weightGraph + " weightAttributes " + weightAttributes + " vertices: " + g.getVertexCount() + " edges: " + g.getEdgeCount());
        numObj = g.getVertexCount();
        numAtt = 1;
        catCount = new CountPack[1];
        attributes = new int[numObj][1];
        int[] counts = new int[singleCat];
        for (int i = 0; i < numObj; i++) {
            attributes[i][0] = singleAtt[i];
            if (singleAtt[i] != -1) {
                counts[singleAtt[i]]++;
            }
        }
        catCount[0] = new CountPack(counts);
        objCoord = new double[numObj][d];
        catCoord = new double[singleCat][d];
        startIndex = new int[1];
        startIndex[0] = 0;
        sumWeights = new double[numObj];
        for (int i = 0; i < numObj; i++) {
            if(! weighted){
            sumWeights[i] += g.getNeighborCount(i) * weightGraph;
            }
            else{
                Collection <MyEdge> e = g.getOutEdges(i);
                for(MyEdge ee: e)
                    sumWeights[i] += ee.getWeight() * weightGraph;
            }
            if (singleAtt[i] != -1) {
                sumWeights[i] += (1 * weightAttributes); //for the categorical attribute
            }
        }

    }

    //init object coordinates randomly and category coords at the center of all assigned objects
    public void init(int seed) {
        Random r = new Random(seed);
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < d; j++) {
                objCoord[i][j] = r.nextDouble();
                for (int k = 0; k < numAtt; k++) {
                    if (attributes[i][k] != -1) {
                        catCoord[startIndex[k] + attributes[i][k]][j] += objCoord[i][j] / (double) catCount[k].counts[attributes[i][k]];
                    }
                }
            }
        }
        //obj();
        displayAtt("init", 0);
        displayGraphAtt(0);
        displayGraphClass();
    }

    private void displayGraphAtt(int attIndex) {
        Visualization v = new Visualization(g);
        int[] a = new int[numObj];
        double[][] cc = new double[numObj][2];
        for (int i = 0; i < numObj; i++) {
            a[i] = attributes[i][attIndex];
            cc[i][0] = objCoord[i][1];
            cc[i][1] = objCoord[i][2];
        }
        if (weighted) {
            v.displayCoordNewWeighted(cc, "Joint Attribute", a);
        } else {
            v.displayCoordNew(cc, "Joint Attribute", a);
        }
    }

    private void displayGraphClass() {
        Visualization v = new Visualization(g);
        double[][] cc = new double[numObj][2];
        for (int i = 0; i < numObj; i++) {
            cc[i][0] = objCoord[i][1];
            cc[i][1] = objCoord[i][2];
        }
        if (weighted) {
            v.displayCoordNewWeighted(cc, "Joint Class", classId);
        } else {
            v.displayCoordNew(cc, "Joint Class", classId);
        }
    }

    private void displayAtt(String title, int attIndex) {
        DataObject[] visu = new DataObject[numObj + catCoord.length];
        for (int i = 0; i < numObj; i++) {
            double[] cc = new double[2];
            for (int j = 0; j < cc.length; j++) {
                cc[j] = objCoord[i][j + 1];
            }
            visu[i] = new DataObject(cc, classId[i], attributes[i][attIndex], i);
        }
        int count = 0;
        for (int i = numObj; i < visu.length; i++) {
            double[] cc = new double[catCoord[count].length - 1];
            for (int j = 0; j < cc.length; j++) {
                cc[j] = catCoord[count][j + 1];
            }
            visu[i] = new DataObject(cc, count + 100, count + 100, i);
            count++;
        }
        VisuAttributes vi = new VisuAttributes(visu, title);
        vi.setSize(600, 600);
        vi.setLocation(200, 0);
        vi.setVisible(true);

//        VisuClass vc = new VisuClass(visu, "class");
//        vc.setSize(600, 600);
//        vc.setLocation(200, 0);
//        vc.setVisible(true);
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
        double attCost = 0.0;
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numAtt; j++) {
                double dist = 0.0;
                for (int l = 0; l < d; l++) {
                    dist += (objCoord[i][l] - catCoord[startIndex[j] + attributes[i][j]][l]) * (objCoord[i][l] - catCoord[startIndex[j] + attributes[i][j]][l]);
                }
                attCost += dist;
            }
        }
        //double cost = weightGraph * graphCost + weightAttributes * attCost;
        //double cost = graphCost;
        double cost = graphCost + attCost;
        System.out.println(graphCost + "  " + attCost + " " + cost);
        return cost;
    }

    //test on synthetic data if this gives the same value
    private double objWeightedNew() {
        double graphCost = 0.0;
        if (weighted) {
            Collection<MyEdge> edges = g.getEdges();
            for (MyEdge e : edges) {
                Pair<Integer> p = g.getEndpoints(e);
                double dist = 0.0;
                for (int l = 0; l < d; l++) {
                    dist += e.getWeight() * (objCoord[p.getFirst()][l] - objCoord[p.getSecond()][l]) * (objCoord[p.getFirst()][l] - objCoord[p.getSecond()][l]);
                }
                graphCost += dist;
            }
        } else {
            Collection<Integer> edges = g.getEdges();
            for (Integer e : edges) {
                Pair<Integer> p = g.getEndpoints(e);
                double dist = 0.0;
                for (int l = 0; l < d; l++) {
                    dist += (objCoord[p.getFirst()][l] - objCoord[p.getSecond()][l]) * (objCoord[p.getFirst()][l] - objCoord[p.getSecond()][l]);
                }
                graphCost += dist;
            }
        }
        double attCost = 0.0;
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numAtt; j++) {
                if (attributes[i][j] != -1) {
                    double dist = 0.0;
                    for (int l = 0; l < d; l++) {
                        dist += (objCoord[i][l] - catCoord[startIndex[j] + attributes[i][j]][l]) * (objCoord[i][l] - catCoord[startIndex[j] + attributes[i][j]][l]);
                    }
                    attCost += dist;
                }
            }
        }
        double cost = weightGraph * graphCost + weightAttributes * attCost;
        System.out.println(weightGraph * graphCost + "  " + weightAttributes * attCost + " " + cost);
        return cost;
    }

    private double objWeighted() {
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

        double attCost = 0.0;
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numAtt; j++) {
                double dist = 0.0;
                for (int l = 0; l < d; l++) {
                    dist += (objCoord[i][l] - catCoord[startIndex[j] + attributes[i][j]][l]) * (objCoord[i][l] - catCoord[startIndex[j] + attributes[i][j]][l]);
                }
                attCost += dist;
            }
        }
        //double cost = weightGraph * graphCost + weightAttributes * attCost;
        //double cost = graphCost;
        double cost = weightGraph * graphCost + weightAttributes * attCost;
        System.out.println(graphCost + "  " + weightAttributes * attCost + " " + cost);
        return cost;
    }

    public void runWeighted() {
        boolean converged = false;
        double minCost = objWeightedNew();
        int iter = 0;
        while (!converged && iter < 100) {
            updateObjectsWeighted();
            updateCategories();
            // displayAtt("first", 0);
            double aktCost = objWeighted();
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
            ea.writeDoubleToMatlab(catCoord, "catCoord", "catCoord");
            double[][] labelsAttribute = new double[numObj][2];
            for (int i = 0; i < numObj; i++) {
                labelsAttribute[i][0] = classId[i];
                labelsAttribute[i][1] = attributes[i][0];
            }
            ea.writeDoubleToMatlab(labelsAttribute, "labelsAttribute", "labelsAttribute");

        }
        displayAtt("finished", 0);
        displayGraphAtt(0);
        displayGraphClass();

    }

    public void run() {
        boolean converged = false;
        double minCost = obj();
        int iter = 0;
        while (!converged && iter < 100) {
            updateObjects();
            updateCategories();
            // displayAtt("first", 0);
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
            ea.writeDoubleToMatlab(catCoord, "catCoord", "catCoord");
            double[][] labelsAttribute = new double[numObj][2];
            for (int i = 0; i < numObj; i++) {
                labelsAttribute[i][0] = classId[i];
                labelsAttribute[i][1] = attributes[i][0];
            }
            ea.writeDoubleToMatlab(labelsAttribute, "labelsAttribute", "labelsAttribute");

        }
        displayAtt("finished", 0);
        displayGraphAtt(0);
        displayGraphClass();

    }

    private void updateCategories() {
        double[][] coordNew = new double[catCoord.length][d];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < d; j++) {
                for (int k = 0; k < numAtt; k++) {
                    coordNew[startIndex[k] + attributes[i][k]][j] += objCoord[i][j] / (double) catCount[k].counts[attributes[i][k]];
                }
            }
        }
        catCoord = coordNew;
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
            for (int k = 0; k < numAtt; k++) {
                for (int j = 0; j < d; j++) {
                    coordNew[i][j] += catCoord[startIndex[k] + attributes[i][k]][j] / sumWeights[i];
                }
            }

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

    private void updateObjectsWeighted() {
        double[][] coordNew = new double[numObj][d];
        for (int i = 0; i < numObj; i++) {
            //graph contribution
            Collection<Integer> n = g.getNeighbors(i);
            for (int p : n) {
                for (int j = 0; j < d; j++) {
                    if (p != i) {
                        coordNew[i][j] += (weightGraph * objCoord[p][j]) / sumWeights[i];
                    }
                }
            }
            //attribute contribution
            for (int k = 0; k < numAtt; k++) {
                for (int j = 0; j < d; j++) {
                    coordNew[i][j] += (weightAttributes * catCoord[startIndex[k] + attributes[i][k]][j]) / sumWeights[i];
                }
            }

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
    
     private void updateObjectsWeightedNew() {
        double[][] coordNew = new double[numObj][d];
        for (int i = 0; i < numObj; i++) {
            //graph contribution
            if(! weighted){
            Collection<Integer> n = g.getNeighbors(i);
            for (int p : n) {
                for (int j = 0; j < d; j++) {
                    if (p != i) {
                        coordNew[i][j] += (weightGraph * objCoord[p][j]) / sumWeights[i];
                    }
                }
            }
            }
            else{
//                Collection<Integer> n = g.getNeighbors(i);
//            for (MyEdge p : n) {
//                for (int j = 0; j < d; j++) {
//                    if (p != i) {
//                        coordNew[i][j] += (weightGraph * objCoord[p][j]) / sumWeights[i];
//                    }
//                }
//            }
                
            }
            //attribute contribution
            for (int k = 0; k < numAtt; k++) {
                for (int j = 0; j < d; j++) {
                    coordNew[i][j] += (weightAttributes * catCoord[startIndex[k] + attributes[i][k]][j]) / sumWeights[i];
                }
            }

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
