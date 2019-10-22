package nature;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import java.util.Collection;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author claudia.plant
 */
public class GraphClustering {

    Graph orig;
    //List<HashMap<Integer, Integer>> ids; //pointNumber -> seedNumber
    HashMap<Integer, Integer>[] pointIdToReachedfrom; //cluster-ids
    HashMap<Integer, Integer>[] pointIdToIndex; //for nodeweights. nodeId -> array index where weight is stored.
    // HashMap<Integer, Integer>[] pointIdToNodename; //Graph: mapping of point id to node name in graph
    Graph[] red;
    double[][] nodeWeights; //numLevels x numNodesPerLevel: how many nodes of the original data are represented by this node
    double[][] edgeWeights; //numLevels x numEdgesPerLevel: how many edgees between the two subsets are in the original data
    int numObj;
    int numLevels;
    int[][] clusterIds;
    int[] numClusters;
    double seedPercentage; //e.g. 0.1: sample 10% of the nodes as seeds
    Random r;

    public GraphClustering(Graph orig, int numLevels, double seedPercentage) {
        this.orig = orig;
        this.numLevels = numLevels;
        this.seedPercentage = seedPercentage;
        // List<HashMap<Integer, Integer>> ids = new ArrayList();
        //HashMap<Integer,Integer>[] ids2 = null ;
        pointIdToReachedfrom = new HashMap[numLevels]; //geht nicht generic Array creation
        pointIdToIndex = new HashMap[numLevels];

        red = new Graph[numLevels];
        numObj = orig.getVertexCount();
        r = new Random(1);
    }

    public double[][] getEdgeWeights() {
        return edgeWeights;
    }

    public double[][] getNodeWeights() {
        return nodeWeights;
    }

    public HashMap<Integer, Integer>[] getPointIdToIndex() {
        return pointIdToIndex;
    }

    public GraphClustering(Graph orig, double seedPercentage) {
        this.orig = orig;
        this.seedPercentage = seedPercentage;
        numLevels = 0;
        double minNumNodesPerLevel = 10.0;
        double expNumNodesLevel0 = ((double) orig.getVertexCount() / 100) * (seedPercentage * 100);
        double expNumNodesCurrentLevel = expNumNodesLevel0;
        double expNumNodesUpperLevel = expNumNodesLevel0;
        while (expNumNodesCurrentLevel > minNumNodesPerLevel) {
            expNumNodesCurrentLevel = expNumNodesUpperLevel / 100 * (seedPercentage * 100);
            numLevels++;
            expNumNodesUpperLevel = expNumNodesCurrentLevel;
        }
        //System.out.println(numLevels);
        // List<HashMap<Integer, Integer>> ids = new ArrayList();
        //HashMap<Integer,Integer>[] ids2 = null ;
        pointIdToReachedfrom = new HashMap[numLevels]; //geht nicht generic Array creation
        pointIdToIndex = new HashMap[numLevels];
//nodeToIndex = new HashMap[numLevels];
        //ids[0] = new HashMap<Double,Integer>() ;
        red = new Graph[numLevels];
        nodeWeights = new double[numLevels][];
        edgeWeights = new double[numLevels][];
        numObj = orig.getVertexCount();
        r = new Random(2); //for uniform, else: 1
        // r = new Random();
        clusterIds = new int[numLevels][orig.getVertexCount()];
        numClusters = new int[numLevels];
    }

    public GraphClustering(Graph orig, double seedPercentage, int seed) {
        this.orig = orig;
        this.seedPercentage = seedPercentage;
        numLevels = 0;
        double minNumNodesPerLevel = 10.0;
        double expNumNodesLevel0 = ((double) orig.getVertexCount() / 100) * (seedPercentage * 100);
        double expNumNodesCurrentLevel = expNumNodesLevel0;
        double expNumNodesUpperLevel = expNumNodesLevel0;
        while (expNumNodesCurrentLevel > minNumNodesPerLevel) {
            expNumNodesCurrentLevel = expNumNodesUpperLevel / 100 * (seedPercentage * 100);
            numLevels++;
            expNumNodesUpperLevel = expNumNodesCurrentLevel;
        }
        //System.out.println(numLevels);
        // List<HashMap<Integer, Integer>> ids = new ArrayList();
        //HashMap<Integer,Integer>[] ids2 = null ;
        pointIdToReachedfrom = new HashMap[numLevels]; //geht nicht generic Array creation
        pointIdToIndex = new HashMap[numLevels];
//nodeToIndex = new HashMap[numLevels];
        //ids[0] = new HashMap<Double,Integer>() ;
        red = new Graph[numLevels];
        nodeWeights = new double[numLevels][];
        edgeWeights = new double[numLevels][];
        numObj = orig.getVertexCount();
        r = new Random(seed); //for uniform, else: 1
        // r = new Random();
    }

    public Graph[] getRed() {
        return red;
    }

    public HashMap<Integer, Integer>[] getIds() {
        return pointIdToReachedfrom;
    }

    public void run() {
        System.out.println("original graph: numVertices: " + orig.getVertexCount() + " numEdges: " + orig.getEdgeCount());
        for (int i = 0; i < numLevels; i++) {
            //createLevel(i);
//            if (i == 9) {
//                System.out.println("m");
//            }
            // createLevelCount(i);
            createLevelCountRandomWalk(i);
            System.out.println("level: " + i + " numVertices: " + red[i].getVertexCount() + " numEdges: " + red[i].getEdgeCount());
        }
        //System.out.println("m");
    }

    //3.3.2017: returns level-wise cluster ids for all vertices. Clusterids start with 0
    public void getClusterIds() {
        for (int i = 0; i < orig.getVertexCount(); i++) {
            clusterIds[0][i] = pointIdToReachedfrom[0].get(i);
            for (int j = 1; j < numLevels; j++) {
                clusterIds[j][i] = pointIdToReachedfrom[j].get(clusterIds[j - 1][i]);
            }
        }
        int[] iid = new int[orig.getVertexCount()];
        for (int i = 0; i < numLevels; i++) {
            for (int j = 0; j < iid.length; j++) {
                iid[j] = -1;
            }
            for (int j = 0; j < orig.getVertexCount(); j++) {
                if (iid[clusterIds[i][j]] == -1) {
                    iid[clusterIds[i][j]] = numClusters[i];
                    numClusters[i]++;
                }
                clusterIds[i][j] = iid[clusterIds[i][j]];
            }
        }
     //   System.out.println("m");

//       for(int i = 0; i < numLevels; i++){
//           //numClusters[i] = pointIdToReachedfrom[i].
//           Object[] idso = pointIdToReachedfrom[i].values().toArray();
//           int[] ids = new int[idso.length];
//           for(int j = 0; j < ids.length; j++){
//               ids[j] = (Integer)idso[j];
//           }
//           HashMap<Integer,Integer> vertexNameToClusterId = new HashMap<Integer,Integer>();
//           for(int j = 0; j < ids.length; j++)
//           vertexNameToClusterId.put(ids[j], j);
//           for(int j = 0; j < orig.getVertexCount(); j++){
//               if(j == 0)
//               System.out.println(j);
//               int id_level_i = pointIdToReachedfrom[i].get(j);
//               int id_upper_level = id_level_i;
//               int id_level = id_level_i;
//               int level = i-1;
//               while(level >=0){
//                   level--;
//                   id_level = pointIdToReachedfrom[level].get(id_upper_level);
//                   id_upper_level = id_level;
//               }
//               clusterIds[i][j] = id_level;
//           }
//       }
    }

    public void createLevelCountRandomWalk(int l) {
        Graph<Integer, Integer> w;
        int lower = l - 1;
        if (lower < 0) {
            w = orig;
        } else {
            w = red[lower];
        }
        Vector<ClusterPath> seeds = new Vector<ClusterPath>();
        pointIdToReachedfrom[l] = new HashMap<Integer, Integer>();
        boolean[] reached = new boolean[w.getVertexCount()];
        Collection<Integer> wv = w.getVertices();
        int countReached = 0;

        HashMap<Integer, Integer> lowerId;
        if (l > 0) {
            lowerId = pointIdToIndex[l - 1];
        } else {
            lowerId = new HashMap<Integer, Integer>();
            int indexCounter = 0;
            for (int currNode : wv) {
                lowerId.put(currNode, indexCounter);
                indexCounter++;
            }
        }

        for (int currNode : wv) {
            double d = r.nextDouble();
            if (d < seedPercentage) {
                ClusterPath cp = new ClusterPath(currNode, currNode);
                seeds.add(cp);
                reached[lowerId.get(currNode)] = true;
                pointIdToReachedfrom[l].put(currNode, currNode);
                countReached++;
            }
        }
        Vector<ClusterPath> activeSeeds = new Vector<ClusterPath>();
        for (int i = 0; i < seeds.size(); i++) {
            activeSeeds.add(seeds.elementAt(i));
        }
        boolean allReached = false;
        boolean restart = false;
        while (!allReached) {
            if (activeSeeds.isEmpty()) {
                for (int i = 0; i < seeds.size(); i++) {
                    activeSeeds.add(seeds.elementAt(i));
                }
                restart = true;
            }
//            if (l == 7 && countReached == 26) {
//                System.out.println("m");
//            }
            int expand = r.nextInt(activeSeeds.size());
            Collection<Integer> n = w.getNeighbors(activeSeeds.elementAt(expand).getPointId());
            Vector<Integer> toExplore = new Vector<Integer>();
            for (int j : n) {
                if (reached[lowerId.get(j)] == false || restart) {
                    toExplore.add(j);
                }
            }
            if (toExplore.isEmpty()) {
                activeSeeds.removeElementAt(expand);
            } else {
                int toFollow = r.nextInt(toExplore.size());
                ClusterPath cp = new ClusterPath(toExplore.elementAt(toFollow), activeSeeds.elementAt(expand).getReachedFrom());
                if (!reached[lowerId.get(toExplore.elementAt(toFollow))]) {
                    pointIdToReachedfrom[l].put(toExplore.elementAt(toFollow), activeSeeds.elementAt(expand).getReachedFrom());
                    reached[lowerId.get(toExplore.elementAt(toFollow))] = true;
                    countReached++;
                }
//                    if(countReached == 135)
//                    System.out.println(countReached);
                activeSeeds.removeElementAt(expand);
                activeSeeds.addElement(cp);

            }
            boolean converged = true;
            for (int i = 0; i < reached.length; i++) {
                if (!reached[i]) {
                    converged = false;
                }
            }
            if (converged) {
                allReached = true;
            }
        }
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();

        HashMap<Integer, Integer> nodeToIndexG = new HashMap<Integer, Integer>();
        int indexCounterG = 0;
        for (int i = 0; i < seeds.size(); i++) {
            g.addVertex(seeds.elementAt(i).getPointId());
            nodeToIndexG.put(seeds.elementAt(i).getPointId(), indexCounterG);
            indexCounterG++;
        }
        double[] gNodeWeights = new double[g.getVertexCount()];

//        for(int i = 0; i < gNodeWeights.length; i++)
//            gNodeWeights[i]  = -1; //substract object itself
        // Collection<Integer> nodeNames = g.getVertices();
        // for (int cn : nodeNames) {
        //count for how many nodes in graph w (level l-1) this node stands
        for (Integer i : pointIdToReachedfrom[l].keySet()) {

            //   if (pointIdToReachedfrom[l].get(i).equals(cn)) {
            //gNodeWeights[nodeToIndexG.get(cn)]++; 
            //add weights of lower levels
//            
//            //TEST//            if(i == 130)
//                System.out.println("m");
//            //TEST
//            
            if (l > 0) {
                gNodeWeights[nodeToIndexG.get(pointIdToReachedfrom[l].get(i))] += nodeWeights[l - 1][pointIdToIndex[l - 1].get(i)];

            } else {
                gNodeWeights[nodeToIndexG.get(pointIdToReachedfrom[l].get(i))]++;
            }
        }
        //}
        //}

        //TEST
        int sumGNodeWeights = 0;
        for (int i = 0; i < gNodeWeights.length; i++) {
            sumGNodeWeights += gNodeWeights[i];
        }
        System.out.println("sumGNodeWeights: " + sumGNodeWeights);
        //TEST

        nodeWeights[l] = gNodeWeights;
        pointIdToIndex[l] = nodeToIndexG;
        //insert all edges as connections in w between all these nodes
        int edgeCounter = 0;
        Vector<Double> edgeMultiplicity = new Vector<Double>();
        for (int currNode : wv) { //all nodes of w (lower level)
            int idCurrent = pointIdToReachedfrom[l].get(currNode);
            Collection<Integer> nBelow = w.getNeighbors(currNode);
            for (int nn : nBelow) {
                int idNeighbor = pointIdToReachedfrom[l].get(nn);
                if (idCurrent < idNeighbor) {
                    if (!g.isNeighbor(idCurrent, idNeighbor)) {
                        g.addEdge(edgeCounter, idCurrent, idNeighbor);
                        int edgeW = w.findEdge(nn, currNode);
                        if (l > 0) {
                            edgeMultiplicity.add(edgeWeights[l - 1][edgeW]);
                        } else {
                            edgeMultiplicity.add(1.0);
                        }
                        edgeCounter++;
                    } else {
                        //increment edgeMultiplicity
                        int edge = g.findEdge(idCurrent, idNeighbor);
                        int edgeW = w.findEdge(nn, currNode);
                        double newMult = 0;
                        if (l > 0) {
                            newMult = edgeMultiplicity.get(edge) + edgeWeights[l - 1][edgeW];
                        } else {
                            newMult = edgeMultiplicity.get(edge) + 1;
                        }
                        edgeMultiplicity.setElementAt(newMult, edge);
                    }
                }
            }
        }
        edgeWeights[l] = new double[edgeMultiplicity.size()];
        for (int i = 0; i < edgeWeights[l].length; i++) {
            edgeWeights[l][i] = edgeMultiplicity.elementAt(i);
        }
        red[l] = g;
    }

    //25.12.15: DEBUG: check if node weight of each node of this level corresponds to the number of nodes in the original graph which each node represents
    private void checkNodeWeights(int l) {
        Graph r = red[l];

    }

    //add hashmap with new ids at the end of the list
    //add reduced Graph of level l: red[l]
    private void createLevel(int l) {
        Graph w;
        int lower = l - 1;
        if (lower < 0) {
            w = orig;
        } else {
            w = red[lower];
        }
        int numVertices = w.getVertexCount();
        Vector<ClusterPath> seeds = new Vector<ClusterPath>();
        pointIdToReachedfrom[l] = new HashMap<Integer, Integer>();
//        pointIdToNodename[l] = new HashMap<Integer, Integer>();
        int counterReached = 0;
        boolean[] reached = new boolean[w.getVertexCount()];

        Collection<Integer> wv = w.getVertices();

        //create index hash map
        HashMap<Integer, Integer> nodeToIndexl = new HashMap<Integer, Integer>();
        int indexCounter = 0;
        for (int currNode : wv) {
            nodeToIndexl.put(currNode, indexCounter);
            indexCounter++;
        }

        for (int currNode : wv) {
            // for (int i = 0; i < numVertices; i++) {
            double d = r.nextDouble();

            if (d < seedPercentage) {
                ClusterPath cp = new ClusterPath(currNode, currNode);
                seeds.add(cp);
                counterReached++;
                reached[nodeToIndexl.get(currNode)] = true;
                //reached[currNode] = true;
                pointIdToReachedfrom[l].put(currNode, currNode);
            }
        }
        // }

//        //check if seeds.size() entries reached
//        int reachedCounter = 0;
//        for(int i = 0; i < reached.length; i++)
//            if(reached[i])
//                reachedCounter++;
//        System.out.println(reachedCounter);
        Vector<ClusterPath> origSeeds = seeds;
//        Object[] origSeeds = seeds.toArray();
        boolean allReached = false;
        int iterCount = 0;
        while (!allReached) {
            Vector<ClusterPath> nextSeeds = new Vector<ClusterPath>();
            for (int i = 0; i < seeds.size(); i++) {
                Collection<Integer> n = w.getNeighbors(seeds.elementAt(i).getPointId());
                for (int j : n) {
                    if (!reached[nodeToIndexl.get(j)]) {
                        counterReached++;
                        ClusterPath cp = new ClusterPath(j, seeds.elementAt(i).getReachedFrom());
                        pointIdToReachedfrom[l].put(j, seeds.elementAt(i).getReachedFrom());
//                        int jj = nodeToIndex.get(j);
//                        if(jj != j)
//                            System.out.println("m");
                        //reached[j] = true;
                        reached[nodeToIndexl.get(j)] = true;
                        nextSeeds.add(cp);
                    }
                }
            }
            boolean converged = true;
            for (int i = 0; i < reached.length; i++) {
                if (!reached[i]) {
                    converged = false;
                    iterCount++;
                }
            }
            if (converged == true) {
                allReached = true;
            } else {
                seeds = nextSeeds;
            }

//            if (counterReached < numVertices) {
//                seeds = nextSeeds;
//            } else {
//                allReached = true;
//            }
        }
        //System.out.println(iterCount);
        //return graph with node names 1...n 
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();

        for (int i = 0; i < origSeeds.size(); i++) {
//            g.addVertex(vertexName);
//            pointIdToNodename[l].put(origSeeds.elementAt(i).getPointId(), vertexName);
//            vertexName++;
            g.addVertex(origSeeds.elementAt(i).getPointId());
        }

        //insert all edges as connections in w between all these nodes
        int edgeCounter = 0;
//        Collection<Integer> wv = w.getVertices();
        for (int currNode : wv) {
            int idCurrent = pointIdToReachedfrom[l].get(currNode);

            Collection<Integer> nBelow = w.getNeighbors(currNode);
            //Collection<Integer> nBelow = w.getNeighbors(pointIdToNodename[indexNodenames].get(currNode));

            for (int nn : nBelow) {

                int idNeighbor = pointIdToReachedfrom[l].get(nn);

                if (idCurrent != idNeighbor && !(g.isNeighbor(idCurrent, idNeighbor))) {

                    g.addEdge(edgeCounter, idCurrent, idNeighbor);

                    edgeCounter++;
                }
            }
        }

        red[l] = g;
    }

    //add hashmap with new ids at the end of the list
    //add reduced Graph of level l: red[l]
    private void createLevelCount(int l) {
        Graph<Integer, Integer> w;
        int lower = l - 1;
        if (lower < 0) {
            w = orig;
        } else {
            w = red[lower];
        }
        Vector<ClusterPath> seeds = new Vector<ClusterPath>();
        pointIdToReachedfrom[l] = new HashMap<Integer, Integer>();
        boolean[] reached = new boolean[w.getVertexCount()];
        Collection<Integer> wv = w.getVertices();
        //create index hash map
        HashMap<Integer, Integer> nodeToIndexW = new HashMap<Integer, Integer>();
        int indexCounter = 0;
        for (int currNode : wv) {
            nodeToIndexW.put(currNode, indexCounter);
            indexCounter++;
        }
        for (int currNode : wv) {
            double d = r.nextDouble();
            if (d < seedPercentage) {
                ClusterPath cp = new ClusterPath(currNode, currNode);
                seeds.add(cp);
                reached[nodeToIndexW.get(currNode)] = true;
                pointIdToReachedfrom[l].put(currNode, currNode);
            }
        }

        //reorder seeds
        int sSize = seeds.size();
        Vector<ClusterPath> origSeeds = new Vector<ClusterPath>();
        Random r = new Random(1);
        boolean[] used = new boolean[sSize];
        for (int i = 0; i < sSize; i++) {
            int index = r.nextInt(sSize);
            while (used[index]) {
                index = r.nextInt(sSize);
            }
            origSeeds.add(seeds.elementAt(index));
            used[index] = true;
        }
        //reorder seeds

        origSeeds = seeds;
        boolean allReached = false;
        // boolean[] reached = new boolean[w.getVertexCount()]; 
        while (!allReached) {
            Vector<ClusterPath> nextSeeds = new Vector<ClusterPath>();
            for (int i = 0; i < seeds.size(); i++) {
                Collection<Integer> n = w.getNeighbors(seeds.elementAt(i).getPointId());
//                 Object[] nn =  n.toArray();
//                 int sCounter = 0;
//                  if (nn.length == 0)
//                    System.out.println("m");
//                int index = r.nextInt(nn.length);
//               
//                sCounter++;
//                while(reached[nodeToIndexW.get((Integer)nn[index])] && sCounter < nn.length){
//                      index = r.nextInt(nn.length);
//                      sCounter++;
//                }

                for (int j : n) {

                    if (!reached[nodeToIndexW.get(j)]) {
                        //  counterReached++;
                        ClusterPath cp = new ClusterPath(j, seeds.elementAt(i).getReachedFrom());
                        pointIdToReachedfrom[l].put(j, seeds.elementAt(i).getReachedFrom());
//                        int jj = nodeToIndex.get(j);
//                        if(jj != j)
//                            System.out.println("m");
                        //reached[j] = true;
                        reached[nodeToIndexW.get(j)] = true;
                        nextSeeds.add(cp);

//                 
                        //  }
                        // }
                    }
                }
            }
            boolean converged = true;
            for (int i = 0; i < reached.length; i++) {
                if (!reached[i]) {
                    converged = false;
                }
            }
            if (converged == true) {
                allReached = true;
            } else {
                seeds = nextSeeds;

                //reorder seeds
                seeds = new Vector<ClusterPath>();
                sSize = nextSeeds.size();
                used = new boolean[sSize];
                for (int i = 0; i < sSize; i++) {
                    int index = r.nextInt(sSize);
                    while (used[index]) {
                        index = r.nextInt(sSize);
                    }
                    seeds.add(nextSeeds.elementAt(index));
                    used[index] = true;
                }
//reorder seeds
            }
        }
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();

        HashMap<Integer, Integer> nodeToIndexG = new HashMap<Integer, Integer>();
        int indexCounterG = 0;
        for (int i = 0; i < origSeeds.size(); i++) {
            g.addVertex(origSeeds.elementAt(i).getPointId());
            nodeToIndexG.put(origSeeds.elementAt(i).getPointId(), indexCounterG);
            indexCounterG++;
        }
        double[] gNodeWeights = new double[g.getVertexCount()];

//        for(int i = 0; i < gNodeWeights.length; i++)
//            gNodeWeights[i]  = -1; //substract object itself
        // Collection<Integer> nodeNames = g.getVertices();
        // for (int cn : nodeNames) {
        //count for how many nodes in graph w (level l-1) this node stands
        for (Integer i : pointIdToReachedfrom[l].keySet()) {
            //   if (pointIdToReachedfrom[l].get(i).equals(cn)) {
            //gNodeWeights[nodeToIndexG.get(cn)]++; 
            //add weights of lower levels
            if (l > 0) {
                gNodeWeights[nodeToIndexG.get(pointIdToReachedfrom[l].get(i))] += nodeWeights[l - 1][nodeToIndexW.get(i)];
            } else {
                gNodeWeights[nodeToIndexG.get(pointIdToReachedfrom[l].get(i))]++;
            }
        }
        //}
        //}

        //TEST
        int sumGNodeWeights = 0;
        for (int i = 0; i < gNodeWeights.length; i++) {
            sumGNodeWeights += gNodeWeights[i];
        }
        System.out.println("sumGNodeWeights: " + sumGNodeWeights);
        //TEST

        nodeWeights[l] = gNodeWeights;
        //insert all edges as connections in w between all these nodes
        int edgeCounter = 0;
        Vector<Double> edgeMultiplicity = new Vector<Double>();
        for (int currNode : wv) { //all nodes of w (lower level)
            int idCurrent = pointIdToReachedfrom[l].get(currNode);
            Collection<Integer> nBelow = w.getNeighbors(currNode);
            for (int nn : nBelow) {
                int idNeighbor = pointIdToReachedfrom[l].get(nn);
                if (idCurrent < idNeighbor) {
                    if (!g.isNeighbor(idCurrent, idNeighbor)) {
                        g.addEdge(edgeCounter, idCurrent, idNeighbor);
                        int edgeW = w.findEdge(nn, currNode);
                        if (l > 0) {
                            edgeMultiplicity.add(edgeWeights[l - 1][edgeW]);
                        } else {
                            edgeMultiplicity.add(1.);
                        }
                        edgeCounter++;
                    } else {
                        //increment edgeMultiplicity
                        int edge = g.findEdge(idCurrent, idNeighbor);
                        int edgeW = w.findEdge(nn, currNode);
                        double newMult = 0;
                        if (l > 0) {
                            newMult = edgeMultiplicity.get(edge) + edgeWeights[l - 1][edgeW];
                        } else {
                            newMult = edgeMultiplicity.get(edge) + 1;
                        }
                        edgeMultiplicity.setElementAt(newMult, edge);
                    }
                }
            }
        }
        edgeWeights[l] = new double[edgeMultiplicity.size()];
        for (int i = 0; i < edgeWeights[l].length; i++) {
            edgeWeights[l][i] = edgeMultiplicity.elementAt(i);
        }
        red[l] = g;
    }
//      //add hashmap with new ids at the end of the list
//    //add reduced Graph of level l: red[l]
//    //counts links and not links between supernodes of level l-1
//    private void createLevelLinkCount(int l) {
//        Graph w;
//        int lower = l - 1;
//        if (lower < 0) {
//            w = orig;
//        } else {
//            w = red[lower];
//        }
//        int numVertices = w.getVertexCount();
//        Vector<ClusterPath> seeds = new Vector<ClusterPath>();
//        pointIdToReachedfrom[l] = new HashMap<Integer, Integer>();
////        pointIdToNodename[l] = new HashMap<Integer, Integer>();
//        int counterReached = 0;
//        boolean[] reached = new boolean[w.getVertexCount()];
//
//        Collection<SuperNode> wv = w.getVertices();
//
//        //create index hash map
//        HashMap<Integer, Integer> nodeToIndex = new HashMap<Integer, Integer>();
//        int indexCounter = 0;
//        for (SuperNode a : wv) {
//            nodeToIndex.put(a.id, indexCounter);
//            indexCounter++;
//        }
//
//
//        for (SuperNode currNode : wv) {
//            // for (int i = 0; i < numVertices; i++) {
//            double d = r.nextDouble();
//
//            if (d < seedPercentage) {
//                ClusterPath cp = new ClusterPath(currNode.id, currNode.id);
//                seeds.add(cp);
//                counterReached++;
//                reached[nodeToIndex.get(currNode.id)] = true;
//                //reached[currNode] = true;
//                pointIdToReachedfrom[l].put(currNode.id, currNode.id);
//            }
//        }
//        // }
//
////        //check if seeds.size() entries reached
////        int reachedCounter = 0;
////        for(int i = 0; i < reached.length; i++)
////            if(reached[i])
////                reachedCounter++;
////        System.out.println(reachedCounter);
//
//        Vector<ClusterPath> origSeeds = seeds;
////        Object[] origSeeds = seeds.toArray();
//        boolean allReached = false;
//        int iterCount = 0;
//        while (!allReached) {
//            Vector<ClusterPath> nextSeeds = new Vector<ClusterPath>();
//            for (int i = 0; i < seeds.size(); i++) {
//                Collection<SuperNode> n = w.getNeighbors(seeds.elementAt(i).getPointId());
//                for (int j : n) {
//                    if (!reached[nodeToIndex.get(j)]) {
//                        counterReached++;
//                        ClusterPath cp = new ClusterPath(j, seeds.elementAt(i).getReachedFrom());
//                        pointIdToReachedfrom[l].put(j, seeds.elementAt(i).getReachedFrom());
////                        int jj = nodeToIndex.get(j);
////                        if(jj != j)
////                            System.out.println("m");
//                        //reached[j] = true;
//                        reached[nodeToIndex.get(j)] = true;
//                        nextSeeds.add(cp);
//                    }
//                }
//            }
//            boolean converged = true;
//            for (int i = 0; i < reached.length; i++) {
//                if (!reached[i]) {
//                    converged = false;
//                    iterCount++;
//                }
//            }
//            if (converged == true) {
//                allReached = true;
//            } else {
//                seeds = nextSeeds;
//            }
//
////            if (counterReached < numVertices) {
////                seeds = nextSeeds;
////            } else {
////                allReached = true;
////            }
//
//        }
//        //System.out.println(iterCount);
//        //return graph with node names 1...n 
//        Graph<SuperNode, SuperEdge> g = new UndirectedSparseGraph<SuperNode, SuperEdge>();
//
//        for (int i = 0; i < origSeeds.size(); i++) {
////            g.addVertex(vertexName);
////            pointIdToNodename[l].put(origSeeds.elementAt(i).getPointId(), vertexName);
////            vertexName++;
//            g.addVertex(new SuperNode(origSeeds.elementAt(i).getPointId(), 0));
//        }
//
//        //insert all edges as connections in w between all these nodes
//        int edgeCounter = 0;
////        Collection<Integer> wv = w.getVertices();
//        for (SuperNode currNode : wv) {
//            int idCurrent = pointIdToReachedfrom[l].get(currNode.id);
//            
//
//            Collection<Integer> nBelow = w.getNeighbors(currNode); //neighbors of currNode at lower level, number of representatives
//            //Collection<Integer> nBelow = w.getNeighbors(pointIdToNodename[indexNodenames].get(currNode));
//
//
//            for (int nn : nBelow) {
//
//
//                int idNeighbor = pointIdToReachedfrom[l].get(nn);
//
//
//
//                if (idCurrent != idNeighbor && !(g.isNeighbor(currNode, idNeighbor))) {
//
//                    g.addEdge(edgeCounter, idCurrent, idNeighbor);
//
//
//                    edgeCounter++;
//                }
//            }
//        }
//      
//        red[l] = g;
//    }
//    
}
