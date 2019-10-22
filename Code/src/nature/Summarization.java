/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import java.util.Collection;
import java.util.Iterator;
import java.util.Vector;

/**
 *
 * @author claudia
 */
public class Summarization {

    Graph g;
    int numObj;
    double[][] dist;
    Vector<SNode> snodes;
    double bestMdl;
    int[] repId;
    double[] viDist; //vi dist of each node to its representant
    double[] entropy; //entropy of all single nodes
    boolean debug = true;

    public Summarization(Graph g) {
        this.g = g;
        numObj = g.getVertexCount();
        dist = new double[numObj][numObj];
        snodes = new Vector<SNode>();
        entropy = new double[numObj];
        this.repId = new int[numObj]; //representant id
        this.viDist = new double[numObj]; // dist to closest representant
        for (int i = 0; i < numObj; i++) //singleton nodes
        {
            repId[i] = i;
            viDist[i] = 0.0;
        }
    }

    public void graphToDistanceMatrix() {
        for (int i = 0; i < dist.length; i++) {
            for (int j = 0; j < dist.length; j++) {
                dist[i][j] = -1.0;
            }
        }
        for (int i = 0; i < dist.length; i++) {
            for (int j = 0; j < dist.length; j++) {
                if (i == j) {
                    dist[i][j] = 0.0;
                } else {
                    if (dist[i][j] == -1.0) {
                        dist[i][j] = vi(i, j);
                        dist[j][i] = dist[i][j];
                    }
                }
            }
        }
        if (debug) {
            IO ea = new IO();
            ea.writeDoubleToMatlab(dist, "dist");
        }

    }

    private double mdl(Vector<SNode> nodes) {
        double mdl = 0.0;
        double idCost = 0.0;
        double codingCost = 0.0;
        for (int i = 0; i < nodes.size(); i++) {
            codingCost += nodes.elementAt(i).entropy;
            codingCost += nodes.elementAt(i).vi;
            idCost += lg2(numObj / nodes.elementAt(i).numAss);
        }
        mdl = idCost + codingCost;
        if (debug) {
            System.out.println(nodes.size() + " nodes, mdl: " + mdl);
        }
        return mdl;
    }

    public void greedyMerge() {
        bestMdl = mdl(snodes);
        boolean improvement = true;
        while (improvement) {
            //search for the index of the next pair of nodes to merge
            int m_i = -1;
            int m_j = -1;
            double minDist = Double.MAX_VALUE;
            for (int i = 0; i < numObj; i++) {
                for (int j = i + 1; j < numObj; j++) {
                    if (dist[i][j] < minDist) {
                        minDist = dist[i][j];
                        m_i = i;
                        m_j = j;
                    }
                }
            }
            boolean bla = merge(m_i, m_j);
            System.out.println("m");

        }

    }

    private boolean merge(int i, int j) {
        //tentatively merge SNodes and check for improvement
        SNode merged = new SNode(snodes.elementAt(i), snodes.elementAt(j));
        //tentatively relable repId; compute for all assigned objects vi
        int removedId = Math.max(i, j);
        int newId = Math.min(i, j);
        int[] repId_n = new int[repId.length];
        for (int ii = 0; ii < numObj; ii++) {
            repId_n[ii] = -1;
        }
        for (int ii = 0; ii < repId_n.length; ii++) {
            if (repId[ii] == removedId) {
                repId_n[ii] = newId;
            } else {
                repId_n[ii] = repId[ii];
            }
        }
        double[] viDist_n = new double[numObj];
        double sumVi = 0.0;
        for (int ii = 0; ii < numObj; ii++) {
            if (repId_n[ii] == newId) {
                viDist_n[ii] = vi(ii, merged);
                sumVi += viDist_n[ii];
            }
        }
        merged.vi = sumVi;
        //compute new mdl
        SNode[] m = new SNode[snodes.size()];
        snodes.copyInto(m);

        for (int ii = 0; ii < m.length; ii++) {
            if (m[ii].label == removedId || m[ii].label == newId) {
                m[ii] = null;
            }
        }
        Vector<SNode> nn = new Vector<SNode>();
        //construct new vector and insert merged at place newId
        for (int ii = 0; ii < m.length; ii++) {
            if (m[ii] != null) {
                nn.addElement(m[ii]);
            } else {
                if (ii == newId) {
                    nn.addElement(merged);
                }
            }


        }
        double newMdl = mdl(nn);
        if (newMdl < bestMdl) {
            bestMdl = newMdl;
            snodes = nn;
            repId = repId_n;
            viDist = viDist_n;
            return true;
        } else {
            return false;
        }
    }

    private double lg2(double d) {
        return Math.log(d) / Math.log(2);
    }

    //vi of a single node to a SNode
    private double vi(int i, SNode s) {
        //compute joint probability of links and remotes required for mi: s.neighbors and remotes if i has link, zero otherwise
        double jpn = 0.0;
        double jpr = 0.0;
        double pLink_s = 0.0;
        double pNotLink_s = 0.0;
        for (int ii = 0; ii < numObj; ii++) {
            pLink_s += s.neighbors[ii];
            pNotLink_s += s.remotes[ii];
            if (g.isNeighbor(i, ii)) {
                jpn += s.neighbors[ii];
            } else {
                jpr += s.remotes[ii];
            }
        }
        jpn++;
        int neighbors_i = g.getNeighborCount(i) + 1; //self link
        int notLinks_i = numObj - neighbors_i;
        double pLink_i = (double) neighbors_i / (double) numObj;
        double pNotLink_i = (double) notLinks_i / (double) numObj;

        pLink_s /= (double) numObj;
        pNotLink_s /= (double) numObj;

        double pCommonLink = jpn / (double) numObj;
        double pCommonNotLink = jpr / (double) numObj;
        double mi = 0.0;
        if (pCommonLink > 0 && pCommonNotLink > 0) {
            mi = pCommonLink * lg2(pCommonLink / (pLink_i * pLink_s)) + pCommonNotLink * lg2(pCommonNotLink / (pNotLink_i * pNotLink_s));
        }
        if (pCommonLink > 0 && pCommonNotLink == 0) {
            mi = pCommonLink * lg2(pCommonLink / (pLink_i * pLink_s));
        }
        if (pCommonNotLink > 0 && pCommonLink == 0) {
            mi = pCommonNotLink * lg2(pCommonNotLink / (pNotLink_i * pNotLink_s));
        }
        double vi = entropy[i] + s.entropy - 2 * mi;
        if (Double.isNaN(vi)) {
            System.err.println("m");
        }
        //System.out.println(vi);
        return vi;


    }

    //vi of single nodes. Also calls the constructor of SNode for singleton nodes i and j
    private double vi(int i, int j) {
        //entropy of i and j
        int neighbors_i = g.getNeighborCount(i) + 1; //self link
        int notLinks_i = numObj - neighbors_i;
        int neighbors_j = g.getNeighborCount(j) + 1; //self link
        int notLinks_j = numObj - neighbors_j;
        double pLink_i = (double) neighbors_i / (double) numObj;
        double pLink_j = (double) neighbors_j / (double) numObj;
        double pNotLink_i = (double) notLinks_i / (double) numObj;
        double pNotLink_j = (double) notLinks_j / (double) numObj;

        double entropy_i = -(pLink_i * lg2(pLink_i) + pNotLink_i * lg2(pNotLink_i));
        double entropy_j = -(pLink_j * lg2(pLink_j) + pNotLink_j * lg2(pNotLink_j));

        int commonNeighbors = 0;
        int commonNotLinks = 0;

        double[] pNeighbors_i = new double[numObj];
        double[] pNeighbors_j = new double[numObj];
        double[] pRemotes_i = new double[numObj];
        double[] pRemotes_j = new double[numObj];
        pNeighbors_i[i] = 1.0;
        pNeighbors_j[j] = 1.0; //self links

        for (int ii = 0; ii < numObj; ii++) {
            if (ii != i && ii != j) {
                boolean neighbor_i = g.isNeighbor(i, ii);
                boolean neighbor_j = g.isNeighbor(j, ii);
                if (neighbor_i) {
                    pNeighbors_i[ii] = 1.0;
                } else {
                    pRemotes_i[ii] = 1.0;
                }
                if (neighbor_j) {
                    pNeighbors_j[ii] = 1.0;
                } else {
                    pRemotes_j[ii] = 1.0;
                }
                if (neighbor_i && neighbor_j) {
                    commonNeighbors++;
                }
                if (!neighbor_i && !neighbor_j) {
                    commonNotLinks++;
                }
            }
        }
        commonNeighbors = commonNeighbors + 2; //self links
        double pCommonLink = (double) commonNeighbors / (double) numObj;
        double pCommonNotLink = (double) commonNotLinks / (double) numObj;
        double mi = 0.0;
        if (pCommonLink > 0 && pCommonNotLink > 0) {
            mi = pCommonLink * lg2(pCommonLink / (pLink_i * pLink_j)) + pCommonNotLink * lg2(pCommonNotLink / (pNotLink_i * pNotLink_j));
        }
        if (pCommonLink > 0 && pCommonNotLink == 0) {
            mi = pCommonLink * lg2(pCommonLink / (pLink_i * pLink_j));
        }
        if (pCommonNotLink > 0 && pCommonLink == 0) {
            mi = pCommonNotLink * lg2(pCommonNotLink / (pNotLink_i * pNotLink_j));
        }
        boolean contains_i = false;
        boolean contains_j = false;
        for (int ii = 0; ii < snodes.size(); ii++) {
            if (snodes.elementAt(ii).label == i) {
                contains_i = true;
            }
            if (snodes.elementAt(ii).label == j) {
                contains_j = true;
            }
        }
        if (!contains_i) {
            SNode si = new SNode(i, pNeighbors_i, pRemotes_i, 1, entropy_i, 0.0);
            entropy[i] = entropy_i;
            snodes.add(si);
        }
        if (!contains_j) {
            SNode sj = new SNode(j, pNeighbors_j, pRemotes_j, 1, entropy_j, 0.0);
            entropy[j] = entropy_j;
            snodes.add(sj);
        }

        double vi = entropy_i + entropy_j - 2 * mi;
        if (Double.isNaN(vi)) {
            System.err.println("m");
        }
        //System.out.println(vi);
        return vi;
    }
}
