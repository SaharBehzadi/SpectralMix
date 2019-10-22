/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.Vector;

/**
 *
 * @author plantc59cs
 */
public class IncrementalEmbedding {

    Graph g;
    Graph e;
    boolean[] embedded;
    Vector<Integer> added; //nodes to be added in the next iteration
    Vector<Integer> active;
    HashMap<Integer, Integer> GtoE; //Hashmaps for node names
    HashMap<Integer, Integer> EtoG;
    //int[] nodeNames; //at index i: name of original node i in graph e
    int currNode; // name of the next node in e
    int currEdge; //name of the next edge in e
    double[][] coord; //numObj x 2
    double[][] result;
    static double convConstant = 1E-3;
    boolean verbose = true;

    public IncrementalEmbedding(Graph g) {
        this.g = g;
        active = new Vector<Integer>();
        added = new Vector<Integer>();
        embedded = new boolean[g.getVertexCount()];
        GtoE = new HashMap<Integer, Integer>();
        EtoG = new HashMap<Integer, Integer>();
        e = new UndirectedSparseGraph<Integer, Integer>();
        currNode = 0;
        currEdge = 0;
        result = new double[g.getVertexCount()][2];
    }

    public void run() {
        Random r = new Random(1);
        //starting vertex: vertex with the largest degree
        int[] degree = new int[g.getVertexCount()];
        int maxDegree = -1;
        int maxIndex = -1;
        for (int i = 0; i < g.getVertexCount(); i++) {
            degree[i] = g.inDegree(i);
            if (degree[i] > maxDegree) {
                maxDegree = degree[i];
                maxIndex = i;
            }
        }
        //add startvertex to e
        e.addVertex(currNode);
        GtoE.put(maxIndex, currNode);
        EtoG.put(currNode, maxIndex);
        currNode++;
        active.add(maxIndex);
        expand(maxIndex);
        coord = new double[e.getVertexCount()][2];
        for (int i = 0; i < e.getVertexCount(); i++) {
            coord[i][0] = r.nextDouble();
            coord[i][1] = r.nextDouble();
        }
        double maxS = 1.0;
//        WeightedMajorization wm = new WeightedMajorization(e, 2, maxS, coord);
//        wm.run(convConstant);
//        double[][] coordOld = wm.coord;

        Grid gg = new Grid(e, maxS);
        gg.run(coord, convConstant);
        double[][] coordOld = gg.getCoord();
        if (verbose) {
            Visualization v = new Visualization(e);
            String s = "result";
            //v.displayCoordNew(coordOld, s);
        }

        while (active.size() > 0 && (e.getVertexCount() < g.getVertexCount())) {
            boolean nodesAdded = false;
            int index = r.nextInt(active.size());
            nodesAdded = expand(active.elementAt(index));
            while (!nodesAdded) {
                index = r.nextInt(active.size());
                nodesAdded = expand(active.elementAt(index));
            }
            //new nodes have been added
            coord = new double[e.getVertexCount()][2];
            for (int i = 0; i < coordOld.length; i++) {
                coord[i][0] = coordOld[i][0];
                coord[i][1] = coordOld[i][1];
            }
            for (int i = 0; i < added.size(); i++) {
                coord[added.elementAt(i)][0] = coordOld[index][0];
                coord[added.elementAt(i)][1] = coordOld[index][1];
            }
//            //DEBUG
//            if (e.getVertexCount() == 151) {
//                IO ea = new IO();
//                String s = "init151";
//                ea.writeDoubleToMatlab(coord, s, "init151");
//            }
//
//            //DEBUG
            gg = new Grid(e, maxS);
            gg.run(coord, convConstant);
            coordOld = gg.getCoord();
            result = gg.getCoord();


//            wm = new WeightedMajorization(e, 2, maxS, coord);
//            wm.run(convConstant);
//            coordOld = wm.coord;
//            result = wm.coord;

            if (verbose) {
                Visualization v = new Visualization(e);
                String s = "result";
                //v.displayCoordNew(coordOld, s);
            }
            // }

            System.out.println(e.getVertexCount() + " " + gg.sigma + " " + gg.getBestCost());
//                System.out.println(e.getVertexCount() + " " + wm.sigma + " " + wm.bestCost);
//            if (e.getVertexCount() == 149 || e.getVertexCount() == 151 || e.getVertexCount() == 138) {
//                IO ea = new IO();
//                String s = "i_" + new Integer(e.getVertexCount()).toString();
//                ea.writeDoubleToMatlab(result, s);
//                String gs = "e_" + new Integer(e.getVertexCount()).toString();
//                ea.writeGraphToMatlab(e, gs);
//            }

        }
         if (verbose) {
                Visualization v = new Visualization(e);
                String s = "result";
                v.displayCoordNew(result, s);
            }
        IO ea = new IO();
        ea.writeDoubleToMatlab(result, "result");

    }

    private void expandNew(int nodenameG) {
        added.removeAllElements();
        Vector<Integer> toAdd = new Vector<Integer>();
        toAdd.add(nodenameG);
        Collection<Integer> nn = g.getNeighbors(nodenameG);
        for (Iterator iterator = nn.iterator(); iterator.hasNext();) {
            Integer index = (Integer) iterator.next();
            toAdd.add(index);
        }
        e = new UndirectedSparseGraph<Integer, Integer>();
        //add nodes
        int nc = 0;
        int ec = 0;
        for (int i = 0; i < toAdd.size(); i++) {
            e.addVertex(nc);
            GtoE.put(toAdd.elementAt(i), nc);
            EtoG.put(nc, toAdd.elementAt(i));
            nc++;
        }
        //add edges between these nodes
        for (int i = 0; i < toAdd.size(); i++) {
            for (int j = 0; j < toAdd.size(); j++) {
                if (g.isNeighbor(EtoG.get(i), EtoG.get(j))) {
                    e.addEdge(ec, i, j);
                }

            }
        }
    }

    //add all neighbors of nodenameG to the graph E; insert all edges between already existing nodes and the new nodes.
    //returns true if new nodes have been inserted
    private boolean expand(int nodenameG) {
        added.removeAllElements();
        Collection<Integer> nn = g.getNeighbors(nodenameG);
        for (Iterator iterator = nn.iterator(); iterator.hasNext();) {
            Integer index = (Integer) iterator.next(); //neighbors of nodenameG in G
            if (!GtoE.containsKey(index)) { //this node needs to be inserted to E
                e.addVertex(currNode);
                GtoE.put(index, currNode);
                EtoG.put(currNode, index);
                added.add(currNode);
                e.addEdge(currEdge, GtoE.get(nodenameG), currNode);
                currEdge++;
                currNode++;
                active.add(index);
            }
        }
        if (added.isEmpty()) {
            return false;
        } else {
            //add all edges between the new nodes and the nodes that are already in e
            for (int i = 0; i < added.size(); i++) {
                Collection<Integer> nnG = g.getNeighbors(EtoG.get(added.elementAt(i)));
                for (Iterator iterator = nnG.iterator(); iterator.hasNext();) {
                    int next = (Integer) iterator.next();
                    if (GtoE.containsKey(next) && !(added.contains(GtoE.get(next)))) { //node has been in e already in the last step
                        if (!(e.isNeighbor(added.elementAt(i), GtoE.get(next)))) {
                            e.addEdge(currEdge, added.elementAt(i), GtoE.get(next));
                            currEdge++;
                        }
                    }

                }
            }
        }
        return true;
    }

}
