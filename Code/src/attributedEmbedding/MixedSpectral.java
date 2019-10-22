/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package attributedEmbedding;

import static attributedEmbedding.JointEmbedding.weighted;
import static attributedEmbedding.JointEmbedding.writeResult;
import datautils.MyEdge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import java.util.Collection;
import java.util.Random;
import java.util.Vector;
import nature.IO;
import nature.Visualization;

/**
 *
 * @author plantc59cs
 */
public class MixedSpectral {

    Graph<Integer, MyEdge>[] g;
    int[][] attributes; //n x numAtt; each categorical attribute value encoded as integer starting from 0
    int numAtt; //number of categorical attributes 
    int numCat; // number of categories
    CountPack[] catCount; //for each categorical attribute the number of objects belonging to each category
    double[][] catCoord; //numCat x d: coords of category centers
    int d; //embedding dimensionality
    double[][] objCoord; //numVertices x d
    int[] classId;
    int[] startIndex; //numAtt: starting index of each Attribute in catCoord
    double[] sumWeights; //for each object: the number of neighbors plus the number of categories; missing values possible
    int numObj;
    double[] weightFactors; //numGraphs + numAtt 
    boolean[] weighted; //numGraphs
    boolean[] directed;
    static boolean writeResult = true;
    static boolean verbose = true;

    public MixedSpectral(Graph[] g, boolean[] weighted, int[][] attributes, int numAtt, int d, int[] classId) {
        this.g = g;
        this.attributes = attributes;
        this.numAtt = numAtt;
        this.d = d;
        this.classId = classId;
        this.weighted = weighted;
        catCount = new CountPack[numAtt];
        numObj = g[0].getVertexCount();
        numCat = 0;
        for (int i = 0; i < numAtt; i++) {
            Vector<Integer> cc = new Vector<Integer>();
            for (int j = 0; j < numObj; j++) {
                if (!cc.contains(attributes[j][i])) {
                    cc.add(attributes[j][i]);
                }
            }
            int numCati = cc.size();
            numCat += numCati;
            int[] counts = new int[numCati];
            for (int j = 0; j < numObj; j++) {
                if (attributes[j][i] != -1) {
                    counts[attributes[j][i]]++;
                }
            }
            catCount[i] = new CountPack(counts);
        }
        objCoord = new double[numObj][d];
        catCoord = new double[numCat][d];
        startIndex = new int[numAtt];
        int placesBefore = 0;
        for (int i = 1; i < startIndex.length; i++) {
            placesBefore += catCount[i - 1].counts.length;
            startIndex[i] = placesBefore;
        }
        //set weight factors for all modalities: all have equal weights. The modality with maximum weight has factor 1, all others are upweighted accordingly
        weightFactors = new double[g.length + numAtt];
        double[] overallWeight = new double[g.length + numAtt];
        double maxWeight = 0.0;
        int maxIndex = -1;
        for (int i = 0; i < g.length; i++) {
            Collection<MyEdge> e = g[i].getEdges();
            for (MyEdge ee : e) {
                overallWeight[i] += ee.getWeight();
            }
            if (overallWeight[i] > maxWeight) {
                maxWeight = overallWeight[i];
                maxIndex = i;
            }
        }
        for (int i = 0; i < numAtt; i++) {
            for (int j = 0; j < catCount[i].counts.length; j++) {
                overallWeight[g.length + i] += catCount[i].counts[j];
            }
            if (overallWeight[g.length + i] > maxWeight) {
                maxWeight = overallWeight[g.length + i];
                maxIndex = g.length + i;
            }
        }
        for (int i = 0; i < weightFactors.length; i++) {
            weightFactors[i] = overallWeight[maxIndex] / overallWeight[i];
        }
        sumWeights = new double[numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < g.length; j++) {
                Collection c = g[j].getSuccessors(i);
                sumWeights[i] += (c.size() * weightFactors[j]); //check this for directed graphs
                //sumWeights[i] += c.size();  //check this for directed graphs
            }
            for (int j = 0; j < numAtt; j++) {
                if (attributes[i][j] != -1) {
                    sumWeights[i] += (1 * weightFactors[g.length + j]); //for the categorical attribute
                    // sumWeights[i] ++; //for the categorical attribute
                }
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
        double d = obj();
        if (verbose) {
            System.out.println("cost after init: " + d);
            displayAtt("sex", 0);
            
            displayGraphClass(0);
            //displayGraphClass(1);
        }
    }

    public void run() {
        boolean converged = false;
        double minCost = obj();
        int iter = 0;
        while (!converged || iter < 100) {
            updateObjects();
            updateCategories();
            // displayAtt("first", 0);
            double aktCost = obj();
            System.out.println(iter + " " + aktCost);
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
        }
        if (verbose) {
            displayAtt("sex", 0);
            displayAtt("schoolclass", 1);
            displayGraphClass(0);
            displayClass();
           // displayGraphClass(1);
        }
    }

    private void updateCategories() {
        double[][] coordNew = new double[catCoord.length][d];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < d; j++) {
                for (int k = 0; k < numAtt; k++) {
                    if (attributes[i][k] != -1) {
                        coordNew[startIndex[k] + attributes[i][k]][j] += objCoord[i][j] / (double) catCount[k].counts[attributes[i][k]];
                    }
                }
            }
        }
        catCoord = coordNew;
    }

    private void updateObjects() {
        double[][] coordNew = new double[numObj][d];
        for (int i = 0; i < numObj; i++) {
            //graph contribution
            for (int j = 0; j < g.length; j++) {
                Collection<Integer> n = g[j].getSuccessors(i);
                for (int p : n) {
                    MyEdge e = g[j].findEdge(i, p);
                    for (int l = 0; l < d; l++) {
                        if (p != i) {
                            coordNew[i][l] += (weightFactors[j] * e.getWeight() * objCoord[p][l]) / sumWeights[i];
                        }
                    }
                }
            }
            //attribute contribution
            for (int k = 0; k < numAtt; k++) {
                if (attributes[i][k] != -1) {
                    for (int j = 0; j < d; j++) {
                        coordNew[i][j] += (weightFactors[g.length + k] * catCoord[startIndex[k] + attributes[i][k]][j]) / sumWeights[i];
                    }
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
    
     private void displayClass() {
        DataObject[] visu = new DataObject[numObj];
        for (int i = 0; i < numObj; i++) {
            double[] cc = new double[d];
            for (int j = 0; j < cc.length; j++) {
                cc[j] = objCoord[i][j];
            }
            visu[i] = new DataObject(cc, classId[i], classId[i], i);
        }
//        int count = 0;
//        for (int i = numObj; i < visu.length; i++) {
//            double[] cc = new double[catCoord[count].length - 1];
//            for (int j = 0; j < cc.length; j++) {
//                cc[j] = catCoord[count][j + 1];
//            }
//            visu[i] = new DataObject(cc, count + 100, count + 100, i);
//            count++;
//        }
        VisuAttributes vi = new VisuAttributes(visu, "class");
        vi.setSize(600, 600);
        vi.setLocation(200, 0);
        vi.setVisible(true);

//        VisuClass vc = new VisuClass(visu, "class");
//        vc.setSize(600, 600);
//        vc.setLocation(200, 0);
//        vc.setVisible(true);
    }

    private void displayAtt(String title, int attIndex) {
        DataObject[] visu = new DataObject[numObj];
        for (int i = 0; i < numObj; i++) {
            double[] cc = new double[d];
            for (int j = 0; j < cc.length; j++) {
                cc[j] = objCoord[i][j];
            }
            visu[i] = new DataObject(cc, classId[i], attributes[i][attIndex], i);
        }
//        int count = 0;
//        for (int i = numObj; i < visu.length; i++) {
//            double[] cc = new double[catCoord[count].length - 1];
//            for (int j = 0; j < cc.length; j++) {
//                cc[j] = catCoord[count][j + 1];
//            }
//            visu[i] = new DataObject(cc, count + 100, count + 100, i);
//            count++;
//        }
        VisuAttributes vi = new VisuAttributes(visu, title);
        vi.setSize(600, 600);
        vi.setLocation(200, 0);
        vi.setVisible(true);

//        VisuClass vc = new VisuClass(visu, "class");
//        vc.setSize(600, 600);
//        vc.setLocation(200, 0);
//        vc.setVisible(true);
    }

    private void displayGraphAtt(int graphIndex, int attIndex) {
        Visualization v = new Visualization(g[graphIndex]);
        int[] a = new int[numObj];
        double[][] cc = new double[numObj][2];
        for (int i = 0; i < numObj; i++) {
            a[i] = attributes[i][attIndex];
            cc[i][0] = objCoord[i][1];
            cc[i][1] = objCoord[i][2];
        }
        if (weighted[graphIndex]) {
            v.displayCoordNewWeighted(cc, "Joint Attribute", a);
        } else {
            v.displayCoordNew(cc, "Joint Attribute", a);
        }
    }

    private void displayGraphClass(int index) {
        Visualization v = new Visualization(g[index]);
        double[][] cc = new double[numObj][2];
        for (int i = 0; i < numObj; i++) {
            cc[i][0] = objCoord[i][1];
            cc[i][1] = objCoord[i][2];
        }
        if (weighted[index]) {
            v.displayCoordNewWeighted(cc, "Joint Class", classId);
        } else {
            v.displayCoordNew(cc, "Joint Class", classId);
        }
    }

    //test on synthetic data if this gives the same value
    private double obj() {
        double[] cost = new double[g.length + numAtt];
        for (int i = 0; i < g.length; i++) {
            Collection<MyEdge> edges = g[i].getEdges();
            for (MyEdge e : edges) {
                Pair<Integer> p = g[i].getEndpoints(e);
                double dist = 0.0;
                for (int l = 0; l < d; l++) {
                    dist += weightFactors[i] * (e.getWeight() * (objCoord[p.getFirst()][l] - objCoord[p.getSecond()][l]) * (objCoord[p.getFirst()][l] - objCoord[p.getSecond()][l]));
                }
                cost[i] += dist;
            }
        }
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numAtt; j++) {
                if (attributes[i][j] != -1) {
                    double dist = 0.0;
                    for (int l = 0; l < d; l++) {
                        dist += weightFactors[g.length + j] * ((objCoord[i][l] - catCoord[startIndex[j] + attributes[i][j]][l]) * (objCoord[i][l] - catCoord[startIndex[j] + attributes[i][j]][l]));
                    }
                    cost[g.length + j] += dist;
                }
            }
        }
        double sumCost = 0.0;
        for (int i = 0; i < cost.length; i++) {
            sumCost += cost[i];
        }
        return sumCost;
    }

}
