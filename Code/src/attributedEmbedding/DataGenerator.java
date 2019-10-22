/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package attributedEmbedding;

import datautils.MyEdge;
import nature.*;
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.filters.FilterUtils;
import edu.uci.ics.jung.algorithms.generators.random.KleinbergSmallWorldGenerator;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import org.apache.commons.collections15.Factory;

/**
 *
 * @author claudia
 */
public class DataGenerator {

    int alpha; //Kleinberg clustering exponent
    int numVertices;
    Graph g;
    int[] attribute;
    int[] classID;

    public DataGenerator() {
    }

    public DataGenerator(int kExponent, int numVertices) {
        this.alpha = kExponent;
        this.numVertices = numVertices;
    }
    Factory<UndirectedSparseGraph<Integer, Integer>> graphFactory
            = new Factory<UndirectedSparseGraph<Integer, Integer>>() {

        public UndirectedSparseGraph<Integer, Integer> create() {
            return new UndirectedSparseGraph<Integer, Integer>();
        }
    };
    Factory<Integer> edgeFactory = new Factory<Integer>() {

        int i = 0;

        public Integer create() {
            return i++;
        }
    };
    Factory<Integer> vertexFactory = new Factory<Integer>() {

        int i = 0;

        public Integer create() {
            return i++;
        }
    };

    public Graph<Integer, Integer> generateKleinberg() {
        KleinbergSmallWorldGenerator ks = new KleinbergSmallWorldGenerator(graphFactory, vertexFactory, edgeFactory, numVertices, numVertices, alpha);
        Graph<Integer, Integer> res = ks.create();
        return res;

    }

    public Graph<Integer, Integer> generateEmptyGraph(int numNodes) {
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();
        for (int i = 0; i < numNodes; i++) {
            g.addVertex(i);
        }
        return g;
    }

    public void clusteredGraphRandomAttributes() {
        int numNodes = 1000;
        classID = new int[numNodes];
        attribute = new int[numNodes];
        Random r = new Random(1);
        int half = Math.floorDiv(numNodes, 2);
        for (int i = 0; i < numNodes; i++) {
            attribute[i] = r.nextInt(2);
            if (i > half) {
                classID[i] = 1;

            }
        }
//        double agreementAttributeWithClass = 0.9;
//        int differentAttribute = (int) Math.round((1.0 - agreementAttributeWithClass) * (double) numNodes);
//        for (int i = 0; i < differentAttribute; i++) {
//            int index = r.nextInt(numNodes);
//            boolean changed = false;
//            if (attribute[index] == 0) {
//                attribute[index] = 1;
//                changed = true;
//            }
//            if (attribute[index] == 1 & !changed) {
//                attribute[index] = 0;
//            }
//        }
        double edgeProbability = 0.01;
        int maxNumEdges = (numNodes * (numNodes - 1)) / 2;
        double clusterConnectionProb = 0.9;
        int numEdges = (int) Math.round(edgeProbability * maxNumEdges);
        g = new UndirectedSparseGraph<Integer, Integer>();
        for (int i = 0; i < numNodes; i++) {
            g.addVertex(i);
        }
        for (int i = 0; i < numEdges; i++) {
            if (r.nextDouble() < clusterConnectionProb) {
                int node1 = r.nextInt(half);
                int node2 = r.nextInt(half);
                while (node1 == node2) {
                    node2 = r.nextInt(half);
                }
                if (r.nextDouble() < 0.5) {
                    g.addEdge(i, node1, node2);
                } else {
                    g.addEdge(i, node1 + half, node2 + half);
                }
            } else {
                int node1 = r.nextInt(numNodes);
                int node2 = r.nextInt(numNodes);
                while (classID[node1] == classID[node2]) {
                    node2 = r.nextInt(numNodes);
                }
                g.addEdge(i, node1, node2);
            }
        }
        // System.out.println("m");

    }

    public void randomGraphRandomAttributes() {
        int numNodes = 1000;
        classID = new int[numNodes];
        attribute = new int[numNodes];
        Random r = new Random(1);
        int half = Math.floorDiv(numNodes, 2);
        for (int i = 0; i < numNodes; i++) {
            attribute[i] = r.nextInt(2);
            if (i > half) {
                classID[i] = 1;

            }
        }
//        double agreementAttributeWithClass = 0.9;
//        int differentAttribute = (int) Math.round((1.0 - agreementAttributeWithClass) * (double) numNodes);
//        for (int i = 0; i < differentAttribute; i++) {
//            int index = r.nextInt(numNodes);
//            boolean changed = false;
//            if (attribute[index] == 0) {
//                attribute[index] = 1;
//                changed = true;
//            }
//            if (attribute[index] == 1 & !changed) {
//                attribute[index] = 0;
//            }
//        }
        double edgeProbability = 0.1;
        int maxNumEdges = (numNodes * (numNodes - 1)) / 2;
        double clusterConnectionProb = 0.9;
        int numEdges = (int) Math.round(edgeProbability * maxNumEdges);
        g = new UndirectedSparseGraph<Integer, Integer>();
        for (int i = 0; i < numNodes; i++) {
            g.addVertex(i);
        }
        for (int i = 0; i < numEdges; i++) {

            int node1 = r.nextInt(numNodes);
            int node2 = r.nextInt(numNodes);

            g.addEdge(i, node1, node2);
        }

        // System.out.println("m");
    }

    public void disagreementBetweenGraphAndAttributes() {
        int numNodes = 1000;
        classID = new int[numNodes];
        attribute = new int[numNodes];
        Random r = new Random(1);
        int half = Math.floorDiv(numNodes, 2);
        for (int i = 0; i < numNodes; i++) {
            if (i < 250) {
                attribute[i] = 1;
            }
            if (i > half) {
                classID[i] = 1;
                //attribute[i] = 1;
            }
        }
//        double agreementAttributeWithClass = 0.9;
//        int differentAttribute = (int) Math.round((1.0 - agreementAttributeWithClass) * (double) numNodes);
//        for (int i = 0; i < differentAttribute; i++) {
//            int index = r.nextInt(numNodes);
//            boolean changed = false;
//            if (attribute[index] == 0) {
//                attribute[index] = 1;
//                changed = true;
//            }
//            if (attribute[index] == 1 & !changed) {
//                attribute[index] = 0;
//            }
//        }
        double edgeProbability = 0.1;
        int maxNumEdges = (numNodes * (numNodes - 1)) / 2;
        double clusterConnectionProb = 0.9;
        int numEdges = (int) Math.round(edgeProbability * maxNumEdges);
        g = new UndirectedSparseGraph<Integer, Integer>();
        for (int i = 0; i < numNodes; i++) {
            g.addVertex(i);
        }
        for (int i = 0; i < numEdges; i++) {
            if (r.nextDouble() < clusterConnectionProb) {
                int node1 = r.nextInt(half);
                int node2 = r.nextInt(half);
                while (node1 == node2) {
                    node2 = r.nextInt(half);
                }
                if (r.nextDouble() < 0.5) {
                    g.addEdge(i, node1, node2);
                } else {
                    g.addEdge(i, node1 + half, node2 + half);
                }
            } else {
                int node1 = r.nextInt(numNodes);
                int node2 = r.nextInt(numNodes);
                while (classID[node1] == classID[node2]) {
                    node2 = r.nextInt(numNodes);
                }
                g.addEdge(i, node1, node2);
            }
        }
        // System.out.println("m");

    }

    public void agreementBetweenGraphAndAttributes() {
        int numNodes = 1000;
        classID = new int[numNodes];
        attribute = new int[numNodes];
        Random r = new Random(1);
        int half = Math.floorDiv(numNodes, 2);
        for (int i = 0; i < numNodes; i++) {
            if (i > half) {
                classID[i] = 1;
                attribute[i] = 1;
            }
        }
        double agreementAttributeWithClass = 1.0;
        int differentAttribute = (int) Math.round((1.0 - agreementAttributeWithClass) * (double) numNodes);
        for (int i = 0; i < differentAttribute; i++) {
            int index = r.nextInt(numNodes);
            boolean changed = false;
            if (attribute[index] == 0) {
                attribute[index] = 1;
                changed = true;
            }
            if (attribute[index] == 1 & !changed) {
                attribute[index] = 0;
            }
        }
        double edgeProbability = 0.01;
        int maxNumEdges = (numNodes * (numNodes - 1)) / 2;
        double clusterConnectionProb = 0.9;
        int numEdges = (int) Math.round(edgeProbability * maxNumEdges);
        g = new UndirectedSparseGraph<Integer, MyEdge>();
        for (int i = 0; i < numNodes; i++) {
            g.addVertex(i);
        }
        for (int i = 0; i < numEdges; i++) {
            if (r.nextDouble() < clusterConnectionProb) {
                int node1 = r.nextInt(half);
                int node2 = r.nextInt(half);
                while (node1 == node2) {
                    node2 = r.nextInt(half);
                }
                if (r.nextDouble() < 0.5) {
                    g.addEdge(new MyEdge(i, 1.0), node1, node2);
                } else {
                    g.addEdge(new MyEdge(i, 1.0), node1 + half, node2 + half);
                }
            } else {
                int node1 = r.nextInt(numNodes);
                int node2 = r.nextInt(numNodes);
                while (classID[node1] == classID[node2]) {
                    node2 = r.nextInt(numNodes);
                }
                g.addEdge(new MyEdge(i, 1.0), node1, node2);
            }
        }

       // g = connect(g, numEdges);

    }

//each node has degree 3: connect to nearest neighbors and insert some additional connections between clusters manually
    public Graph<Integer, Integer> generateClusteredGraph(int[] ids, double[][] coords) {
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();
        for (int i = 0; i < ids.length; i++) {
            g.addVertex(i);
        }
        int edgeCounter = 0;

        for (int i = 0; i < g.getVertexCount(); i++) {
            double[] dd = new double[g.getVertexCount()];
            dd[i] = Double.MAX_VALUE;
            for (int j = 0; j < g.getVertexCount(); j++) {
                if (i != j) {
                    dd[j] = dist(coords[i], coords[j]);
                }
            }
            double min = Double.MAX_VALUE;
            int minIndex = -1;
            for (int j = 0; j < dd.length; j++) {
                if (dd[j] < min) {
                    min = dd[j];
                    minIndex = j;
                }
            }
            g.addEdge(edgeCounter, i, minIndex); //1NN
            edgeCounter++;
            dd[minIndex] = Double.MAX_VALUE;
            min = Double.MAX_VALUE;
            minIndex = -1;
            for (int j = 0; j < dd.length; j++) {
                if (dd[j] < min) {
                    min = dd[j];
                    minIndex = j;
                }
            }
            g.addEdge(edgeCounter, i, minIndex); //2NN
            edgeCounter++;
            dd[minIndex] = Double.MAX_VALUE;
            min = Double.MAX_VALUE;
            minIndex = -1;
            for (int j = 0; j < dd.length; j++) {
                if (dd[j] < min) {
                    min = dd[j];
                    minIndex = j;
                }
            }
            g.addEdge(edgeCounter, i, minIndex); //3NN
            edgeCounter++;
            dd[minIndex] = Double.MAX_VALUE;
            min = Double.MAX_VALUE;
            minIndex = -1;
            for (int j = 0; j < dd.length; j++) {
                if (dd[j] < min) {
                    min = dd[j];
                    minIndex = j;
                }
            }
            g.addEdge(edgeCounter, i, minIndex); //4NN
            edgeCounter++;
            dd[minIndex] = Double.MAX_VALUE;
            min = Double.MAX_VALUE;
            minIndex = -1;
            for (int j = 0; j < dd.length; j++) {
                if (dd[j] < min) {
                    min = dd[j];
                    minIndex = j;
                }
            }
            g.addEdge(edgeCounter, i, minIndex); //5NN
            edgeCounter++;
            dd[minIndex] = Double.MAX_VALUE;
            min = Double.MAX_VALUE;
            minIndex = -1;
            for (int j = 0; j < dd.length; j++) {
                if (dd[j] < min) {
                    min = dd[j];
                    minIndex = j;
                }
            }
            g.addEdge(edgeCounter, i, minIndex); //6NN
            edgeCounter++;
//            dd[minIndex] = Double.MAX_VALUE;
//            min = Double.MAX_VALUE;
//            minIndex = -1;
//            for (int j = 0; j < dd.length; j++) {
//                if (dd[j] < min) {
//                    min = dd[j];
//                    minIndex = j;
//                }
//            }
//            g.addEdge(edgeCounter, i, minIndex); //4NN
//            edgeCounter++;
//            dd[minIndex] = Double.MAX_VALUE;
//            min = Double.MAX_VALUE;
//            minIndex = -1;
//            for (int j = 0; j < dd.length; j++) {
//                if (dd[j] < min) {
//                    min = dd[j];
//                    minIndex = j;
//                }
//            }
//            g.addEdge(edgeCounter, i, minIndex); //4NN
//            edgeCounter++;
//            dd[minIndex] = Double.MAX_VALUE;
//            min = Double.MAX_VALUE;
//            minIndex = -1;
//            for (int j = 0; j < dd.length; j++) {
//                if (dd[j] < min) {
//                    min = dd[j];
//                    minIndex = j;
//                }
//            }
//            g.addEdge(edgeCounter, i, minIndex); //4NN
//            edgeCounter++;
//            dd[minIndex] = Double.MAX_VALUE;
//            min = Double.MAX_VALUE;
//            minIndex = -1;
//            for (int j = 0; j < dd.length; j++) {
//                if (dd[j] < min) {
//                    min = dd[j];
//                    minIndex = j;
//                }
//            }
//            g.addEdge(edgeCounter, i, minIndex); //4NN
//            edgeCounter++;
        }
        //hand added Edges
        g.addEdge(edgeCounter, 32, 5);
        edgeCounter++;
        g.addEdge(edgeCounter, 162, 147);
        edgeCounter++;
        g.addEdge(edgeCounter, 85, 63);
        edgeCounter++;
        g.addEdge(edgeCounter, 32, 5);
        edgeCounter++;
        g.addEdge(edgeCounter, 102, 115);
        edgeCounter++;
        g.addEdge(edgeCounter, 79, 106);
        edgeCounter++;
        g.addEdge(edgeCounter, 183, 187);
        edgeCounter++;
        g.addEdge(edgeCounter, 40, 79);
        edgeCounter++;
        g.addEdge(edgeCounter, 115, 49);
        edgeCounter++;
        //between clusters
        g.addEdge(edgeCounter, 1, 111);
        edgeCounter++;
        g.addEdge(edgeCounter, 8, 111);
        edgeCounter++;
        g.addEdge(edgeCounter, 24, 47);
        edgeCounter++;
        g.addEdge(edgeCounter, 24, 99);
        edgeCounter++;
        g.addEdge(edgeCounter, 119, 166);
        edgeCounter++;
        g.addEdge(edgeCounter, 165, 114);
        edgeCounter++;
        g.addEdge(edgeCounter, 19, 186);
        edgeCounter++;
        g.addEdge(edgeCounter, 31, 165);
        edgeCounter++;
        g.addEdge(edgeCounter, 30, 9);
        edgeCounter++;
        g.addEdge(edgeCounter, 25, 5);
        edgeCounter++;
        g.addEdge(edgeCounter, 26, 5);
        edgeCounter++;
        return g;
    }

    //connect each pair of objects with distance < t, coord: numObj x dim
    public Graph<Integer, Integer> generateThresholdGraph(double t, int numObj) {
        int dim = 2;
        Random r = new Random(2);
        double[][] coord = new double[numObj][dim];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < dim; j++) {
                coord[i][j] = r.nextDouble();
            }
        }
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();
        int edgeCounter = 0;
        for (int i = 0; i < coord.length; i++) {
            g.addVertex(i);
            for (int j = i + 1; j < coord.length; j++) {
                if (dist(coord[i], coord[j]) < t) {
                    g.addEdge(edgeCounter, i, j);
                    edgeCounter++;
                }
            }
        }
//        double[][] pd = pathdist(g);
//        IO ea = new IO();
//        ea.writeDoubleToMatlab(pd, "pathDist");
//        ea.writeDoubleToMatlab(coord, "coord");
        return g;
    }

    //add edges to make the graph connected
    public Graph connect(Graph g, int edgeCounter) {
        Random r = new Random(1);
        boolean connected = false;
        while (!connected) {
            WeakComponentClusterer<Number, Number> clusterer = new WeakComponentClusterer<Number, Number>();
            Set<Set<Number>> clusterset = clusterer.transform(g);
            if (clusterset.size() > 1) {
                for (Set<Number> cluster : clusterset) {
                    if (cluster.size() == 1) {
                        Object[] nodeName = cluster.toArray();
                        Integer name = (Integer) nodeName[0];
                        int other = r.nextInt(g.getVertexCount());
                        while (other == name) {
                            other = r.nextInt(g.getVertexCount());
                        }
                        g.addEdge(edgeCounter, name, other);
                        edgeCounter++;
                    } else {
                        Object[] nodeName = cluster.toArray();
                        int index = r.nextInt(nodeName.length);
                        Integer name = (Integer) nodeName[index];
                        int other = r.nextInt(g.getVertexCount());
                        while (other == name) {
                            other = r.nextInt(g.getVertexCount());
                        }
                        g.addEdge(edgeCounter, name, other);
                        edgeCounter++;
                    }
                    //System.out.println("m");
                }
            } else {
                connected = true;
            }
        }
        // g = getLargestComponent(g);
        return g;
    }

    public Graph getLargestComponent(Graph g) {
        FilterUtils filt = new FilterUtils();
        WeakComponentClusterer<Number, Number> clusterer = new WeakComponentClusterer<Number, Number>();
        Set<Set<Number>> clusterset = clusterer.transform(g);
        Set<Number> largest = Collections.EMPTY_SET;
        for (Set<Number> cluster : clusterset) {
            if (cluster.size() > largest.size()) {
                largest = cluster;
            }
        }
        Graph<Integer, Integer> res = new UndirectedSparseGraph<Integer, Integer>();
        boolean[] isIncluded = new boolean[g.getVertexCount()];
        Object[] nodeNames = largest.toArray();
        int[] index = new int[nodeNames.length];
        //add all nodes
        for (int i = 0; i < index.length; i++) {
            Integer name = (Integer) nodeNames[i];
            isIncluded[name] = true;
            res.addVertex(i);
            index[i] = i;
        }
        for (int i = 0; i < isIncluded.length; i++) {
            if (!isIncluded[i]) {
                System.out.println("missing node: " + i);
            }
        }
        //add all Edges

        for (int i = 0; i < nodeNames.length; i++) {
            for (int j = 0; j < nodeNames.length; j++) {
                if (i > j && g.isNeighbor(nodeNames[i], nodeNames[j])) {
                    res.addEdge(getIndex(i, j, nodeNames.length), i, j);

                }
            }
        }
//        //check if graph is connected
//        for(int i = 0; i < res.getVertexCount(); i++){
//            int bla = res.getNeighborCount(i);
//           // if(bla == 0)
//                System.out.println(i + " " + bla);
//        }

        return res;
    }

    //without diagonal
    public int getIndex(int index1, int index2, int numObj) {
        int i1 = Math.min(index1, index2);
        int i2 = Math.max(index1, index2);
        int count = i1;
        int dec = 1;
        int placesBefore = 0;
        while (count > 0) {
            placesBefore = placesBefore + (numObj - dec);
            dec++;
            count--;
        }
        int offset = i2 - i1 - 1;
        return placesBefore + offset;
    }

    private double[][] pathdist(Graph<Integer, Integer> g) {
        DijkstraShortestPath<Integer, Integer> alg = new DijkstraShortestPath(g);
        int n = g.getVertexCount();
        double[][] dist = new double[n][n];
        double maxDist_mds = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i > j) {
                    List<Integer> l = alg.getPath(i, j);
                    if (l.size() > 0) {
                        dist[i][j] = l.size();
                        dist[j][i] = dist[i][j];
                        if (dist[i][j] > maxDist_mds) {
                            maxDist_mds = dist[i][j];
                        }
                    } else {
                        dist[i][j] = -1;
                        dist[j][i] = -1;
                    }

                }
            }
        }
        for (int i = 0; i < dist.length; i++) {
            for (int j = 0; j < dist.length; j++) {
                if (dist[i][j] == 0.0 && i != j) {
                    //dist[i][j] = maxDist_mds + 1;
                }
            }
        }
        return dist;
    }

    public Graph<Integer, Integer> generateKleinbergSelf() {
        int size = (int) Math.sqrt(numVertices);
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();
        int latDim = 2;
        int[][] lattice = new int[numVertices][latDim];
        int row = 0;
        int col = 0;
        for (int i = 0; i < numVertices; i++) {
            g.addVertex(i);
            lattice[i][0] = row;
            lattice[i][1] = col;
            col++;
            if (col == size) {
                row++;
                col = 0;
            }
        }
        //add edges to all 4 neighbors in lattice
        int edgeCounter = 0;
        for (int i = 0; i < numVertices; i++) {
            row = lattice[i][0];
            col = lattice[i][1];
            int left = i - 1;
            int right = i + 1;
            int top = i - size;
            int bottom = i + size;
            if (col > 0) {
                g.addEdge(edgeCounter, i, left);
                edgeCounter++;
            }
            if (col < size - 1) {
                g.addEdge(edgeCounter, i, right);
                edgeCounter++;
            }
            if (row > 0) {
                g.addEdge(edgeCounter, i, top);
                edgeCounter++;
            }
            if (row < size - 1) {
                g.addEdge(edgeCounter, i, bottom);
                edgeCounter++;
            }

        }
        // Add long range connections
        int longRangeCounter = 0;
        Random r = new Random(1);
        for (int i = 0; i < numVertices; i++) {
            for (int j = i + 1; j < numVertices; j++) {
                if (!g.isNeighbor(i, j)) {
                    double dist = manhDist(lattice[i], lattice[j]);
                    double p = 1.0 / Math.pow(dist, alpha);
                    double factor = 10;
                    p /= factor;
                    double rr = r.nextDouble();
                    if (rr < p) {
                        g.addEdge(edgeCounter, i, j);
                        edgeCounter++;
                        longRangeCounter++;
                    }
                }
            }
        }
        System.out.println("longRange: " + longRangeCounter);

        return g;

    }

    public static double manhDist(int[] x, int[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += Math.abs(x[i] - y[i]);
        }
        return result;
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
