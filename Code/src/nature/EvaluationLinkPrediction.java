/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.filters.FilterUtils;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StreamTokenizer;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

/**
 *
 * @author plantc59cs
 */
public class EvaluationLinkPrediction {

    Graph train;
    int[] startEdge; //start and endpoints of edges to be predicted
    int[] endEdge;
    int[] startNoEdge; //start and endpoints of not-edges
    int[] endNoEdge;
    double[][] coord;
    int d;
    HashMap<Integer, Integer> nodenameToIndexTrain;
    static int numExperiments = 10;
    static int numTrain = 10000;
    static int numTest = 1000;
    static int toPredict = 5000;

   
    public EvaluationLinkPrediction(double[][] coord, int d) {
        this.coord = coord;
        this.d = d;
        nodenameToIndexTrain = new HashMap<Integer, Integer>();
    }

    public EvaluationLinkPrediction() {
        nodenameToIndexTrain = new HashMap<Integer, Integer>();
    }
    
    
    
    
//     public EvaluationLinkPrediction(int d) {
//        this.coord = new double[]
//        this.d = d;
//        nodenameToIndexTrain = new HashMap<Integer, Integer>();
//    }

    
    
    public void edgesToPredict(){
        String file = "stackoverflowTest.txt";
        int numLinks = 5000; //5000 links to predict
        boolean[] valid = new boolean[numLinks];
        try {
            File source = new File(file);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(source)));
            st1.eolIsSignificant(false);
            st1.parseNumbers();
            int counter = 0;
            int numValid = 0;
//            HashMap<Integer, Integer> nodenameToIndex = new HashMap<Integer, Integer>();
//            int nodeCounter = 0;
            while (st1.ttype != StreamTokenizer.TT_EOF && counter < numLinks) {
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex1 = (int) st1.nval;
                //int v1 = vertex1 - 1;
                
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex2 = (int) st1.nval;
               
                if((vertex1 != vertex2) && train.containsVertex(nodenameToIndexTrain.get(vertex1)) && train.containsVertex(nodenameToIndexTrain.get(vertex2))){
                valid[counter] = true;
                counter++;
                numValid++;
                }
            
            }
            System.out.println(numValid);

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
//writes nodenameToIndexTrain, train
    public void stackoverflowTrainGraph() {
        String file = "stackoverflowTrain.txt";
        int numLinks = 501550; //5000 links to predict
        int[] start = new int[numLinks];
        int[] end = new int[numLinks];
        Graph g = new UndirectedSparseGraph<Integer, Integer>();
        Vector<Integer> nodes = new Vector<Integer>();
        //original nodename to index, i.e. new nodename in train graph
        HashMap<Integer,Integer> nodenameToIndex = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> indexToNodename = new HashMap<Integer, Integer>();
       
        try {
            File source = new File(file);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(source)));
            st1.eolIsSignificant(false);
            st1.parseNumbers();
            int counter = 0;
//            HashMap<Integer, Integer> nodenameToIndex = new HashMap<Integer, Integer>();
//            int nodeCounter = 0;
            while (st1.ttype != StreamTokenizer.TT_EOF && counter < numLinks) {
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex1 = (int) st1.nval;
                //int v1 = vertex1 - 1;
                if (!nodes.contains(vertex1)) {
                    nodes.add(vertex1);
//                    nodenameToIndex.put(v1, nodeCounter);
//                    nodeCounter++;
                }
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex2 = (int) st1.nval;
                if (!nodes.contains(vertex2)) {
                    nodes.add(vertex2);
                }
//                if(counter == 501549)
//                    System.out.println("m");
                start[counter] = vertex1;
                end[counter] = vertex2;
                counter++;
            }
            System.out.println(nodes.size());

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        for (int i = 0; i < nodes.size(); i++) {
            int currNode = nodes.elementAt(i);
            nodenameToIndex.put(currNode, i);
            indexToNodename.put(i, currNode);
            g.addVertex(i);
        }
        int edgeCounter = 0;
        for (int i = 0; i < start.length; i++) {
            if (g.containsVertex(nodenameToIndex.get(start[i])) && g.containsVertex(nodenameToIndex.get(end[i])) && (start[i] != end[i]) && !(g.isNeighbor(nodenameToIndex.get(start[i]), nodenameToIndex.get(end[i])))) {
                g.addEdge(edgeCounter, nodenameToIndex.get(start[i]), nodenameToIndex.get(end[i]));
                edgeCounter++;
            }
        }
        System.out.println("trainGraph: " + g.getVertexCount() + " " + g.getEdgeCount());
        FilterUtils filt = new FilterUtils();
        WeakComponentClusterer<Number, Number> clusterer = new WeakComponentClusterer<Number, Number>();
        Set<Set<Number>> clusterset = clusterer.transform(g);
        Set<Number> largest = Collections.EMPTY_SET;
        for (Set<Number> cluster : clusterset) {
            if (cluster.size() > largest.size()) {
                largest = cluster;
            }
        }
        //HashMap<Integer, Integer>nodenameToIndexNew = new HashMap<Integer, Integer>();
        Graph<Integer, Integer> res = new UndirectedSparseGraph<Integer, Integer>();
        boolean[] isIncluded = new boolean[g.getVertexCount()];
        Object[] nodeNames = largest.toArray();
        int[] index = new int[nodeNames.length];
        //add all nodes
        for (int i = 0; i < index.length; i++) {
            Integer name = (Integer) nodeNames[i];
          //  nodenameToIndexNew.put(nodeNa, name)
            isIncluded[name] = true;
            res.addVertex(i);
            nodenameToIndexTrain.put(indexToNodename.get(name), i);
            index[i] = i;
        }

        int counter = 0;

        for (int i = 0; i < isIncluded.length; i++) {
            if (!isIncluded[i]) {
                System.out.println("missing node: " + i);
            } else {
//                gtCoord[counter][0] = gt[counter][0];
//                gtCoord[counter][1] = gt[counter][1];
                counter++;
            }
        }
       
        //add all Edges
        edgeCounter = 0;
        for (int i = 0; i < nodeNames.length; i++) {
            for (int j = 0; j < nodeNames.length; j++) {
                if (i > j && g.isNeighbor(nodeNames[i], nodeNames[j])) {
                    res.addEdge(edgeCounter, i, j);
                    edgeCounter++;
                }
            }
        }
//        //check if graph is connected
//        for(int i = 0; i < res.getVertexCount(); i++){
//            int bla = res.getNeighborCount(i);
//           // if(bla == 0)
//                System.out.println(i + " " + bla);
//        }

//        return res;     
        System.out.println("largest Component: " + res.getVertexCount() + " " + res.getEdgeCount());
        this.train = res;
        DataUtils du = new DataUtils();
        du.writeAdjacencyList(res, "stackoverflowAdj");
    }

    public void linkPredictionExperimentsHadamard(int numLinks) {
        //create the arff data set
        int numVar = d + 1;
        FastVector fvWekaAttributes = new FastVector(numVar);
        for (int i = 0; i < numVar - 1; i++) {
            Attribute a = new Attribute("d" + i);
            fvWekaAttributes.addElement(a);
        }
        FastVector my_nominal_values = new FastVector(2);
        my_nominal_values.addElement("true");
        my_nominal_values.addElement("false");
        Attribute cl = new Attribute("class", my_nominal_values);
        fvWekaAttributes.addElement(cl);
        Instances train = new Instances("Graph", fvWekaAttributes, numTrain);
        Instances test = new Instances("Graph", fvWekaAttributes, numTest);
        train.setClassIndex(d);
        test.setClassIndex(d);
//        for (int i = 0; i < g.getVertexCount(); i++) {
//            Instance inst = new Instance(nodes.numAttributes());
//            for (int j = 0; j < nodes.numAttributes() - 1; j++) {
//                inst.setValue((Attribute) fvWekaAttributes.elementAt(j), coord[i][j]);
//            }
//            inst.setValue((Attribute) fvWekaAttributes.elementAt(fvWekaAttributes.size() - 1), "c" + labels[i]);
//            nodes.add(inst);
//        }

    }

}
