/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.visualization3d.VisualizationViewer;
import java.util.Collection;
import java.util.HashMap;
import java.util.Random;

/**
 *
 *
 */
public class HierarchicalEmbeddingWeightedMoreD {

    Graph g;
    Graph<Integer, Integer>[] red;
    double[][] nodeWeights;
    double[][] edgeMultiplicities;
    HashMap<Integer, Integer>[] pointIdToReachedfrom;
    double[][][] coordL;
    HashMap<Integer, Integer>[] pointIdToNodename; //for each level: mappint of point id to node name, coord-index in the embedding
    HashMap<Integer, Integer>[] nodenameToPointId;
    HashMap<Integer, Integer>[] pointIdToIndex; // point id -> index in array node weights
    int numLevels;
    double pe;
    boolean verbose = true;
    static double convConstant = 1E-3;
    int d;

    //alpha: edgeWeights; beta: nodeWeightes Endpunkt * nodeWeights anderer Endpunkt - alpha
    public HierarchicalEmbeddingWeightedMoreD(Graph g, double seedPercentage, int d) {
        this.d = d;
        this.g = g;
        Random r = new Random();
        int bla = r.nextInt(100);
        System.out.println(bla);
        GraphClustering gc = new GraphClustering(g, seedPercentage, bla);
        gc.run();
        red = gc.getRed();
        nodeWeights = gc.getNodeWeights();
        edgeMultiplicities = gc.getEdgeWeights();
        pointIdToReachedfrom = gc.getIds();
        pointIdToIndex = gc.getPointIdToIndex();

        numLevels = red.length;
        //relabel the graphs for vertices staring from 0
        pointIdToNodename = new HashMap[numLevels]; //pointId : numbers in original graph -> nodeName: numbers of vertices starting from zero
        nodenameToPointId = new HashMap[numLevels];
        double[][] relabeledMult = new double[numLevels][];
        double[][] relabeledNodeWeights = new double[numLevels][];
        for (int i = 0; i < numLevels; i++) {
            relabeledMult[i] = new double[edgeMultiplicities[i].length];
            relabeledNodeWeights[i] = new double[nodeWeights[i].length];
            pointIdToNodename[i] = new HashMap<Integer, Integer>();
            nodenameToPointId[i] = new HashMap<Integer, Integer>();
            int nodename = 0;
            Collection<Integer> wv = red[i].getVertices();
            for (int currNode : wv) {
                pointIdToNodename[i].put(currNode, nodename);
                nodenameToPointId[i].put(nodename, currNode);
                nodename++;
            }
            Graph rel = new UndirectedSparseGraph<Integer, Integer>();
            for (int j = 0; j < pointIdToNodename[i].size(); j++) {
                rel.addVertex(j);
                int nodeNameOrig = nodenameToPointId[i].get(j);
                int indexWeight = pointIdToIndex[i].get(nodeNameOrig);
                relabeledNodeWeights[i][j] = nodeWeights[i][indexWeight];
            }
            //Collection<Integer> ed = red[i].getEdges();
            int edgeCounter = 0;
            for (int currNode : wv) {
                Collection<Integer> nn = red[i].getNeighbors(currNode);
                int indexCurr = pointIdToNodename[i].get(currNode);
                for (int currN : nn) {
                    int indexNn = pointIdToNodename[i].get(currN);
                    if (!rel.isNeighbor(indexCurr, indexNn)) {
                        rel.addEdge(edgeCounter, indexCurr, indexNn);
                        int indexEdgeOrig = red[i].findEdge(currNode, currN);
                        relabeledMult[i][edgeCounter] = edgeMultiplicities[i][indexEdgeOrig];
                        edgeCounter++;
                    }
                }
            }

            red[i] = rel;
            nodeWeights[i] = relabeledNodeWeights[i];
            edgeMultiplicities[i] = relabeledMult[i];
        }

        coordL = new double[numLevels][][];
    }

    public double[][] returnCoord() {
        int level = red.length - 1;
        while (level >= 0) {
            if (verbose) {
                System.out.println("---------Level: " + level);
            }
            double[][] coord = new double[red[level].getVertexCount()][d];
            //Embedding expects Graph with nodes 0..n
            if (level == red.length - 1) { //random init
                Random r = new Random(5);
                // r = new Random();
                for (int i = 0; i < coord.length; i++) {
                    for (int j = 0; j < d; j++) {
                        coord[i][j] = r.nextDouble();
                    }
                }
            } else { //init coords of each object to the coords of the cluster representative of the upper level, z.B. level 2
                for (int i = 0; i < red[level].getVertexCount(); i++) {
                    for (int j = 0; j < d; j++) {
                        //  int bla = pointIdToNodename[level+1].get(35); //need nodenameToPointID here
                        int father = pointIdToReachedfrom[level + 1].get(nodenameToPointId[level].get(i));
                        int fatherIndex = pointIdToNodename[level + 1].get(father);
                        coord[i][j] = coordL[level + 1][fatherIndex][j];
                    }
                }
                coord = blowUpDisplay(coordL[level + 1], coord);

            }
            double maxS = 1.0;
//            //2D////////////////////////////////////////////////////////////////
            //embedd

            //Grid gg = new Grid(red[level], maxS);
//            //Weighted
            HierarchicalWeightedGrid gg = new HierarchicalWeightedGrid(red[level], maxS, nodeWeights[level], edgeMultiplicities[level]);
//            //Weighted

            IO ea = new IO();
//            if (level == 3 || level == 10 || level == 8) {
//                String filename = "graphLevel" + level + ".mat";
//                ea.writeGraphToMatlab(red[level], filename);
//            }
            //gg.run();
            if (red[level].getVertexCount() > 1) {
                gg.run(coord, convConstant);
                coordL[level] = gg.getCoord();
            } else {
                coordL[level] = new double[1][d];
                Random r = new Random(1);
                for (int i = 0; i < d; i++) {
                    coordL[level][0][i] = r.nextDouble();
                }
            }
            if (verbose) {
                Visualization v = new Visualization(red[level]);
                String s = "result" + level;
                //v.displayCoordNew(coordL[level], s);
            }
//           
//              String filename = "coord_"+level+".mat";
//              ea.writeDoubleToMatlab(coordL[level], "coord", filename);

//            //2D//////////////////////////////////////////////////////////////////////////
//            //3D///////////////////////////////////////////////////////////////////////
//            WeightedMajorization wm = new WeightedMajorization(red[level], 3, 1, maxS);
//            wm.setGroundTruth(coord);
//            wm.groundTruthInit = true;
//            wm.run(convConstant);
//            coordL[level] = wm.coord;
//3D/////////////////////////////////////////////////////////////////////////////////////////////////////
            level--;
        }
        double[][] coord = new double[g.getVertexCount()][d];
        for (int i = 0; i < g.getVertexCount(); i++) {
            for (int j = 0; j < d; j++) {
                //  int bla = pointIdToNodename[level+1].get(35); //need nodenameToPointID here
                int father = pointIdToReachedfrom[0].get(i);
                int fatherIndex = pointIdToNodename[0].get(father);
                coord[i][j] = coordL[0][fatherIndex][j];
            }
        }
        coord = blowUpDisplay(coordL[0], coord);

        double maxS = 20.0;
//        ////////////2D////////////////////////////////////////////////////////////////
        Grid gg = new Grid(g, maxS);
        gg.run(coord, convConstant);
        double[][] cRes = gg.getCoord();
        pe = gg.getBestCost();
        if (verbose) {
            Visualization v = new Visualization(g);
            String s = "result";
           // v.displayCoordNew(cRes, s);
        }
//        ///////////////////////////2D////////////////////////////////////////////////////////

        ////////////////3D//////////////////////////////////////////////////////////
//        DisplayEmbedding3d dd = new DisplayEmbedding3d(g, coord);
//        dd.display();
//        double[][] cRes = dd.getCoord();
        //3D///////////////////////////////////////////////////////////////////////////////
        DataUtils du = new DataUtils();
        du.saveAsMatlab(cRes, "coord", "result.mat");
        return cRes;
    }

    public double[][] run() {
        int level = red.length - 1;
        while (level >= 0) {
            if (verbose) {
                System.out.println("---------Level: " + level);
            }
            double[][] coord = new double[red[level].getVertexCount()][d];
            //Embedding expects Graph with nodes 0..n
            if (level == red.length - 1) { //random init
                Random r = new Random(5);
                // r = new Random();
                for (int i = 0; i < coord.length; i++) {
                    for (int j = 0; j < d; j++) {
                        coord[i][j] = r.nextDouble();
                    }
                }
            } else { //init coords of each object to the coords of the cluster representative of the upper level, z.B. level 2
                for (int i = 0; i < red[level].getVertexCount(); i++) {
                    for (int j = 0; j < d; j++) {
                        //  int bla = pointIdToNodename[level+1].get(35); //need nodenameToPointID here
                        int father = pointIdToReachedfrom[level + 1].get(nodenameToPointId[level].get(i));
                        int fatherIndex = pointIdToNodename[level + 1].get(father);
                        coord[i][j] = coordL[level + 1][fatherIndex][j];
                    }
                }
                coord = blowUpDisplay(coordL[level + 1], coord);

            }
            double maxS = 1.0;
//            //2D////////////////////////////////////////////////////////////////
            //embedd

            //Grid gg = new Grid(red[level], maxS);
//            //Weighted
// public WeightedMajorizationWeighted(Graph g, int d, double maxSigma, double[] nodeWeights, double[] edgeMultiplicity, double[][] init)
            WeightedMajorizationWeighted gg = new WeightedMajorizationWeighted(red[level], d, maxS, nodeWeights[level], edgeMultiplicities[level]);
            //wj.run(1e-3);
            //     HierarchicalWeightedGrid gg = new HierarchicalWeightedGrid(red[level], maxS, nodeWeights[level], edgeMultiplicities[level]);
//            //Weighted

            IO ea = new IO();
//            if (level == 3 || level == 10 || level == 8) {
//                String filename = "graphLevel" + level + ".mat";
//                ea.writeGraphToMatlab(red[level], filename);
//            }
            //gg.run();
            if (red[level].getVertexCount() > 1) {
                gg.run(coord, convConstant);
                coordL[level] = gg.getCoord();
            } else {
                coordL[level] = new double[1][d];
                Random r = new Random(1);
                for (int i = 0; i < d; i++) {
                    coordL[level][0][i] = r.nextDouble();
                }
            }
            if (verbose) {
                Visualization v = new Visualization(red[level]);
                String s = "result" + level;
                //v.displayCoordNew(coordL[level], s);
            }
//           
//              String filename = "coord_"+level+".mat";
//              ea.writeDoubleToMatlab(coordL[level], "coord", filename);

//            //2D//////////////////////////////////////////////////////////////////////////
//            //3D///////////////////////////////////////////////////////////////////////
//            WeightedMajorization wm = new WeightedMajorization(red[level], 3, 1, maxS);
//            wm.setGroundTruth(coord);
//            wm.groundTruthInit = true;
//            wm.run(convConstant);
//            coordL[level] = wm.coord;
//3D/////////////////////////////////////////////////////////////////////////////////////////////////////
            level--;
        }
        double[][] coord = new double[g.getVertexCount()][d];
        for (int i = 0; i < g.getVertexCount(); i++) {
            for (int j = 0; j < d; j++) {
                //  int bla = pointIdToNodename[level+1].get(35); //need nodenameToPointID here
                int father = pointIdToReachedfrom[0].get(i);
                int fatherIndex = pointIdToNodename[0].get(father);
                coord[i][j] = coordL[0][fatherIndex][j];
            }
        }
        coord = blowUpDisplay(coordL[0], coord);

        double maxS = 20.0;
//        ////////////2D////////////////////////////////////////////////////////////////
        WeightedMajorization gg = new WeightedMajorization(g, d, maxS, coord);
        gg.run(convConstant);
        double[][] cRes = gg.coord;
        pe = gg.bestCost;
        if (verbose) {
            Visualization v = new Visualization(g);
            String s = "result";
           // v.displayCoordNew(cRes, s);
        }
//        ///////////////////////////2D////////////////////////////////////////////////////////

        ////////////////3D//////////////////////////////////////////////////////////
//        DisplayEmbedding3d dd = new DisplayEmbedding3d(g, coord);
//        dd.display();
//        double[][] cRes = dd.getCoord();
        //3D///////////////////////////////////////////////////////////////////////////////
        DataUtils du = new DataUtils();
        //du.saveAsMatlab(cRes, "coord", "result.mat");
        return cRes;
    }

    public void run(double tau) {
        int level = red.length - 1;
        while (level >= 0) {
            if (verbose) {
                System.out.println("---------Level: " + level);
            }
            double[][] coord = new double[red[level].getVertexCount()][d];
            //Embedding expects Graph with nodes 0..n
            if (level == red.length - 1) { //random init
                Random r = new Random(5);
                // r = new Random();
                for (int i = 0; i < coord.length; i++) {
                    for (int j = 0; j < d; j++) {
                        coord[i][j] = r.nextDouble();
                    }
                }
            } else { //init coords of each object to the coords of the cluster representative of the upper level, z.B. level 2
                for (int i = 0; i < red[level].getVertexCount(); i++) {
                    for (int j = 0; j < d; j++) {
                        //  int bla = pointIdToNodename[level+1].get(35); //need nodenameToPointID here
                        int father = pointIdToReachedfrom[level + 1].get(nodenameToPointId[level].get(i));
                        int fatherIndex = pointIdToNodename[level + 1].get(father);
                        coord[i][j] = coordL[level + 1][fatherIndex][j];
                    }
                }
                coord = blowUpDisplay(coordL[level + 1], coord);

            }
            double maxS = 1.0;
//            //2D////////////////////////////////////////////////////////////////
            //embedd

            Grid gg = new Grid(red[level], maxS);
            IO ea = new IO();
//            if (level == 3 || level == 10 || level == 8) {
//                String filename = "graphLevel" + level + ".mat";
//                ea.writeGraphToMatlab(red[level], filename);
//            }
            //gg.run();
            gg.setTau(tau);
            gg.run(coord, convConstant);
            coordL[level] = gg.getCoord();
            if (verbose) {
                Visualization v = new Visualization(red[level]);
                String s = "result" + level;
                //v.displayCoordNew(coordL[level], s);
            }
//           
//              String filename = "coord_"+level+".mat";
//              ea.writeDoubleToMatlab(coordL[level], "coord", filename);

//            //2D//////////////////////////////////////////////////////////////////////////
//            //3D///////////////////////////////////////////////////////////////////////
//            WeightedMajorization wm = new WeightedMajorization(red[level], 3, 1, maxS);
//            wm.setGroundTruth(coord);
//            wm.groundTruthInit = true;
//            wm.run(convConstant);
//            coordL[level] = wm.coord;
//3D/////////////////////////////////////////////////////////////////////////////////////////////////////
            level--;
        }
        double[][] coord = new double[g.getVertexCount()][d];
        for (int i = 0; i < g.getVertexCount(); i++) {
            for (int j = 0; j < d; j++) {
                //  int bla = pointIdToNodename[level+1].get(35); //need nodenameToPointID here
                int father = pointIdToReachedfrom[0].get(i);
                int fatherIndex = pointIdToNodename[0].get(father);
                coord[i][j] = coordL[0][fatherIndex][j];
            }
        }
        coord = blowUpDisplay(coordL[0], coord);

        double maxS = 20.0;
//        ////////////2D////////////////////////////////////////////////////////////////
        Grid gg = new Grid(g, maxS);
        gg.setTau(tau);
        gg.run(coord, convConstant);
        double[][] cRes = gg.getCoord();
        pe = gg.getBestCost();
        if (verbose) {
            Visualization v = new Visualization(g);
            String s = "result";
            //v.displayCoordNew(cRes, s);
        }
//        ///////////////////////////2D////////////////////////////////////////////////////////

        ////////////////3D//////////////////////////////////////////////////////////
//        DisplayEmbedding3d dd = new DisplayEmbedding3d(g, coord);
//        dd.display();
//        double[][] cRes = dd.getCoord();
        //3D///////////////////////////////////////////////////////////////////////////////
        DataUtils du = new DataUtils();
        du.saveAsMatlab(cRes, "coord", "result.mat");
    }

    public double[][] blowUpDisplay(double[][] coordOld, double[][] coordNew) {
        int numObjOld = coordOld.length;
        int numObjNew = coordNew.length;
        double[][] minMax = new double[d][2]; //0: minimum, 1: maximum
        for (int i = 0; i < d; i++) {
            minMax[i][0] = Double.MAX_VALUE;
            minMax[i][1] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < numObjOld; i++) {
            for (int j = 0; j < d; j++) {
                if (coordOld[i][j] < minMax[j][0]) {
                    minMax[j][0] = coordOld[i][j];
                }
                if (coordOld[i][j] > minMax[j][1]) {
                    minMax[j][1] = coordOld[i][j];
                }
            }
        }
        //double[] spread = new double[d];
        double[] scaleFactor = new double[d];
        for (int j = 0; j < d; j++) {
            double spread = minMax[j][1] - minMax[j][0];
            double spacePerObjOld = spread / (double) numObjOld;
            double spreadNew = spread + (numObjNew - numObjOld) * spacePerObjOld;
            scaleFactor[j] = spreadNew / spread;

        }
        for (int i = 0; i < numObjNew; i++) {
            for (int j = 0; j < d; j++) {

                coordNew[i][j] *= scaleFactor[j];
            }
        }
        return coordNew;

    }
}
