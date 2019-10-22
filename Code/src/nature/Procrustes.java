/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;
import java.util.Collection;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;

/**
 *
 * @author plantc59cs
 */
public class Procrustes {

    Graph g;
    int numSeeds;
    int[] clId;
    int n;
    boolean[] clusterBorder;
    HashMap<Integer, Integer>[] origToSmall; //contains for each cluster the Hashmap OrigNodeNames -> SmallNodeNames
    HashMap<Integer, Integer>[] smallToOrig;
    Vector<Integer>[] sharedNodes;
    Graph[] s; //contains for each cluster the small graphs
    double[][][] coord;
    static double convConstant = 1E-3;
    boolean verbose = true;

    public Procrustes(Graph g, int numSeeds) {
        this.g = g;
        this.numSeeds = numSeeds;
        this.n = g.getVertexCount();
        this.clId = new int[n];
        for (int i = 0; i < n; i++) {
            clId[i] = -1;
        }
        clusterBorder = new boolean[n];
        s = new Graph[numSeeds];
        origToSmall = new HashMap[numSeeds];
        smallToOrig = new HashMap[numSeeds];
        sharedNodes = new Vector[numSeeds * (numSeeds - 1) / 2];
        for (int i = 0; i < sharedNodes.length; i++) {
            sharedNodes[i] = new Vector<Integer>();
        }
        coord = new double[numSeeds][][];
    }

    public void embedd() {
        for (int i = 0; i < numSeeds; i++) {
            double maxS = 20.0;
            Grid gg = new Grid(s[i], maxS);
            gg.run(convConstant);
            coord[i] = gg.getCoord();
            if (verbose) {
                int[] clId_i = new int[s[i].getVertexCount()];
                for (int k = 0; k < clId_i.length; k++) {
                    clId_i[k] = clId[smallToOrig[i].get(k)];
                }
                Visualization v = new Visualization(s[i]);
                v.displayCoordNew(coord[i], new Integer(i).toString(), clId_i);
            }
        }
    }

    public void mergeClusters(int a, int b) {
        //get all shared nodes and perform a Procrustes analysis    
        //write the joint graph to s[Math.min(a,b)]
        //write the joint coord to coord[Math.min(a,b)][][]
        DataUtils du = new DataUtils();
        int clusterPairIndex = du.getIndex(a, b, numSeeds);
        int numShared = sharedNodes[clusterPairIndex].size();

        double[][] sharedCoord_a = new double[numShared][2];
        double[][] sharedCoord_b = new double[numShared][2];
        for (int i = 0; i < numShared; i++) {
            int index_a = origToSmall[a].get(sharedNodes[clusterPairIndex].elementAt(i));
            int index_b = origToSmall[b].get(sharedNodes[clusterPairIndex].elementAt(i));
            sharedCoord_a[i][0] = coord[a][index_a][0];
            sharedCoord_a[i][1] = coord[a][index_a][1];
            sharedCoord_b[i][0] = coord[a][index_b][0];
            sharedCoord_b[i][1] = coord[a][index_b][1];
        }
        //DEBUG
        du.saveAsMatlab2(sharedCoord_a, sharedCoord_b, "sharedCoord_a", "sharedCoord_b", "inputProc.mat");

        //DEBUG
        double[] mean_a = new double[2];
        double[] mean_b = new double[2];
        for (int i = 0; i < numShared; i++) {
            mean_a[0] += sharedCoord_a[i][0] / (double) numShared;
            mean_a[1] += sharedCoord_a[i][1] / (double) numShared;
            mean_b[0] += sharedCoord_b[i][0] / (double) numShared;
            mean_b[1] += sharedCoord_b[i][1] / (double) numShared;
        }
        double var_a = 0.0;
        double var_b = 0.0;
        for (int i = 0; i < numShared; i++) {
            var_a += (sharedCoord_a[i][0] - mean_a[0]) * (sharedCoord_a[i][0] - mean_a[0]);
            sharedCoord_a[i][0] -= mean_a[0];
            var_a += (sharedCoord_a[i][1] - mean_a[1]) * (sharedCoord_a[i][1] - mean_a[1]);
            sharedCoord_a[i][1] -= mean_a[1];
            var_b += (sharedCoord_b[i][0] - mean_b[0]) * (sharedCoord_b[i][0] - mean_b[0]);
            sharedCoord_b[i][0] -= mean_b[0];
            var_b += (sharedCoord_b[i][1] - mean_b[1]) * (sharedCoord_b[i][1] - mean_b[1]);
            sharedCoord_b[i][1] -= mean_b[1];
        }
        var_a /= (numShared * 2);
        var_a = Math.sqrt(var_a);
        var_b /= (numShared * 2);
        var_b = Math.sqrt(var_b);

        for (int i = 0; i < numShared; i++) {
            sharedCoord_a[i][0] /= var_a;
            sharedCoord_a[i][1] /= var_a;
            sharedCoord_b[i][0] /= var_b;
            sharedCoord_b[i][1] /= var_b;
        }
        du.saveAsMatlab2(sharedCoord_a, sharedCoord_b, "sharedCoord_a", "sharedCoord_b", "afterScaling.mat");
       
        double upper = 0.0;
        double lower = 0.0;
        for (int i = 0; i < numShared; i++) {
            upper += sharedCoord_b[i][0] * sharedCoord_a[i][1] - sharedCoord_b[i][1] * sharedCoord_a[i][0];
            lower += sharedCoord_b[i][0] * sharedCoord_a[i][0] - sharedCoord_b[i][1] * sharedCoord_a[i][1];
        }
        double phi = Math.atan(upper / lower);
        
//DEBUG
        double[][] landmarks_a = new double[numShared][2];
        double[][] landmarks_b = new double[numShared][2];
        for (int i = 0; i < numShared; i++) {
            landmarks_a[i][0] = (coord[a][i][0] - mean_a[0]) / var_a;
            landmarks_a[i][1] = (coord[a][i][1] - mean_a[1]) / var_a;
        }
        for (int i = 0; i < numShared; i++) {
            double[] c_nor = new double[2];
            c_nor[0] = (coord[b][i][0] - mean_b[0]) / var_b;
            c_nor[1] = (coord[b][i][1] - mean_b[1]) / var_b;
            landmarks_b[i][0] = Math.cos(phi) * c_nor[0] - Math.sin(phi) * c_nor[1];
            landmarks_b[i][1] = Math.sin(phi) * c_nor[0] + Math.cos(phi) * c_nor[1];
        }

//        du.saveAsMatlab(jointCoord, "jointCoord", "resultOfProc.mat");
        du.saveAsMatlab2(landmarks_a, landmarks_b, "landmarks_a", "landmarks_b", "landmarks.mat");
        
      
        //DEBUG
        //write new coord array, shared nodes contained only once
        //nodes of a plus shared nodes: not modified
        //nodes of b except shared nodes between a and b rotated
        double[][] jointCoord = new double[s[a].getVertexCount() + s[b].getVertexCount() - numShared][2];
        HashMap<Integer, Integer> origToSmall_m = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> smallToOrig_m = new HashMap<Integer, Integer>();
        Graph<Integer, Integer> m = new UndirectedSparseGraph<Integer, Integer>();
        //all from a
        for (int i = 0; i < coord[a].length; i++) {
            jointCoord[i][0] = (coord[a][i][0] - mean_a[0]) / var_a;
            jointCoord[i][1] = (coord[a][i][1] - mean_a[1]) / var_a;
            m.addVertex(i);
            origToSmall_m.put(smallToOrig[a].get(i), i);
            smallToOrig_m.put(i, smallToOrig[a].get(i));
        }
        //all from b except shared nodes and border nodes shared by both
        int counter = coord[a].length;
        for (int i = 0; i < coord[b].length; i++) {
            //shared node
            int origName = smallToOrig[b].get(i);

            if (!sharedNodes[clusterPairIndex].contains(smallToOrig[b].get(i)) && !origToSmall_m.containsKey(origName)) {
                double[] c_nor = new double[2];
                c_nor[0] = (coord[b][i][0] - mean_b[0]) / var_b;
                c_nor[1] = (coord[b][i][1] - mean_b[1]) / var_b;
                jointCoord[counter][0] = Math.cos(phi) * c_nor[0] - Math.sin(phi) * c_nor[1];
                jointCoord[counter][1] = Math.sin(phi) * c_nor[0] + Math.cos(phi) * c_nor[1];
                m.addVertex(counter);
                origToSmall_m.put(smallToOrig[b].get(i), counter);
                smallToOrig_m.put(counter, smallToOrig[b].get(i));
                counter++;
            }
        }
        //add edges to the graph
        //add all edges of a
        int edgeCounter = 0;
        Collection<Integer> edges_a = s[a].getEdges();
        for (int j : edges_a) {
            Pair<Integer> ep = s[a].getEndpoints(j);
            int name_p1Orig = smallToOrig[a].get(ep.getFirst());
            int name_p2Orig = smallToOrig[a].get(ep.getSecond());
            if (name_p1Orig == 125 || name_p2Orig == 125) {
                System.out.println("m");
            }
            //   if(!m.isIncident(origToSmall_m.get(name_p1Orig), origToSmall_m.get(name_p1Orig))){
            m.addEdge(edgeCounter, origToSmall_m.get(name_p1Orig), origToSmall_m.get(name_p2Orig));
            edgeCounter++;
            // }

        }
        //add all edges of b except shared nodes
        Collection<Integer> edges_b = s[b].getEdges();
        for (int j : edges_b) {
            Pair<Integer> ep = s[b].getEndpoints(j);
            int name_p1Orig = smallToOrig[b].get(ep.getFirst());
            int name_p2Orig = smallToOrig[b].get(ep.getSecond());
            if (!m.isIncident(origToSmall_m.get(name_p1Orig), origToSmall_m.get(name_p2Orig))) {
                m.addEdge(edgeCounter, origToSmall_m.get(name_p1Orig), origToSmall_m.get(name_p2Orig));
                edgeCounter++;
            }
        }
        //DEBUG: Find node without edges
        for (int i = 0; i < m.getVertexCount(); i++) {
            Collection<Integer> n = m.getNeighbors(i);
            if (n.isEmpty()) {
                System.out.println(i + " " + i + " origNode: " + smallToOrig_m.get(i));
            }
        }
        Visualization v = new Visualization(m);
        v.getCoordinatesItMajWithInit(jointCoord);
        // v.displayCoordNew(jointCoord, "proc");

        System.out.println("m");
    }

    //TO DO: manage the active seeds in a vector
    public void cluster() {
        Random r = new Random(1);
        Vector<Integer> activeSeeds = new Vector<Integer>();
        boolean finished = false;
        for (int i = 0; i < numSeeds; i++) {
            int seedId = r.nextInt(n);
            while (clId[seedId] != -1) {
                seedId = r.nextInt(n);
            }
            clId[seedId] = i;
            activeSeeds.add(seedId);
        }
        while (!finished) {
            int toExpand = r.nextInt(activeSeeds.size());
            Collection<Integer> n = g.getNeighbors(activeSeeds.elementAt(toExpand));
            int currCl = clId[activeSeeds.elementAt(toExpand)];
            activeSeeds.remove(toExpand);
            for (int j : n) {
                if (clId[j] == -1) {
                    clId[j] = currCl;
                    activeSeeds.add(j);
                }
            }
            if (activeSeeds.isEmpty()) {
                finished = true;
            }

        }

        if (verbose) {
            Visualization v = new Visualization(g);
            v.displayCoordNew(v.getCoordinatesIsomapOnly(), "clusteringIso", clId);
        }

    }

    private void indexTest() {
        DataUtils du = new DataUtils();
        for (int i = 0; i < numSeeds; i++) {
            for (int j = 0; j < numSeeds; j++) {
                if (i != j) {
                    int[] p = du.getPair(du.getIndex(i, j, numSeeds), numSeeds);
                    System.out.println("i " + i + " j " + j + " " + " index " + du.getIndex(i, j, numSeeds) + " pair: " + +p[0] + " " + p[1]);

                    // System.out.println("index: " + p[0] + " " + p[1]);
                }
            }
        }
    }

    //15.3.2017: patition the graph for embedding. Writes s and l. In each small graph all border nodes are included
    public void partition() {
        //indexTest();
        //identify the border nodes
        DataUtils du = new DataUtils();
        for (int i = 0; i < g.getVertexCount(); i++) {
            if (!clusterBorder[i]) {
                Collection<Integer> n = g.getNeighbors(i);
                for (int j : n) {
                    if (clId[j] != clId[i]) {
                        clusterBorder[i] = true;
                        clusterBorder[j] = true;
                    }
                }
            }
        }
        //construct small graphs for each cluster
        for (int i = 0; i < numSeeds; i++) {
            //System.out.println("seed: " + i);
            s[i] = new UndirectedSparseGraph<Integer, Integer>();
            origToSmall[i] = new HashMap<Integer, Integer>();
            smallToOrig[i] = new HashMap<Integer, Integer>();

            //include all vertices of clId i plus the border nodes of i to another cluster
            int nodeCounter = 0;
            for (int j = 0; j < n; j++) {
                if (clId[j] == i) {
                    s[i].addVertex(nodeCounter);
                    origToSmall[i].put(j, nodeCounter);
                    smallToOrig[i].put(nodeCounter, j);
                    nodeCounter++;
                }

                if (clId[j] != i && clusterBorder[j]) {
                    boolean neighborsIn_i = false;
                    Collection<Integer> n = g.getNeighbors(j);
                    for (int k : n) {
                        if (clId[k] == i) {
                            neighborsIn_i = true;
                        }
                    }
                    if (neighborsIn_i) {
                        s[i].addVertex(nodeCounter);
                        if (!sharedNodes[du.getIndex(i, clId[j], numSeeds)].contains(j)) {
                            sharedNodes[du.getIndex(i, clId[j], numSeeds)].add(j);
                        }
                        origToSmall[i].put(j, nodeCounter);
                        smallToOrig[i].put(nodeCounter, j);
                        nodeCounter++;
                    }
                }
            }
            //add all edges between any pair of nodes which are included
            int edgeCounter = 0;
            for (int j = 0; j < s[i].getVertexCount(); j++) {
//                if(j == 64)
//                System.out.println(j);
                //if j is an inner node of the original graph: include all edges
                int origName_j = smallToOrig[i].get(j);
                if (clId[origName_j] == i && !clusterBorder[origName_j]) {
                    Collection<Integer> n = g.getNeighbors(origName_j);
                    for (int k : n) {
//                        if(k == 57)
//                        System.out.println("k " + k);
                        s[i].addEdge(edgeCounter, j, origToSmall[i].get(k));
                        edgeCounter++;
                    }
                }
                //if j is a border node: include only edges to the cluster and to other border nodes
                if (clId[origName_j] != i && clusterBorder[origName_j]) {
                    Collection<Integer> n = g.getNeighbors(origName_j);
                    for (int k : n) {
                        if (smallToOrig[i].containsValue(k)) {
                            s[i].addEdge(edgeCounter, j, origToSmall[i].get(k));
                            edgeCounter++;
                        }
                    }
                }
            }

        }
    }
}
