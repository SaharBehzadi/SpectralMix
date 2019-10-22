/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import java.util.Vector;

/**
 *
 * @author claudia
 */
public class SNode {

    int label; //this is used as key identifier
    //Vector<Integer> objects;
    double[] neighbors; //fractional link pattern of all nodes to this SNode
    double[] remotes;
    int numAss; //number of assigned objects
    double entropy;
    double vi; //vi

    public SNode(int label, double[] neighbors, double[] remotes, int numAss, double entropy, double vi) {
        this.label = label;
        this.numAss = numAss;
        this.neighbors = neighbors;
        this.remotes = remotes;
        this.entropy = entropy;
        this.vi = vi;
    }

    public SNode(SNode i, SNode j) {
        this.label = Math.min(i.label, j.label);
        this.neighbors = new double[i.neighbors.length];
        this.remotes = new double[i.remotes.length];
        this.numAss = i.numAss + j.numAss;
        double n_n = 0.0; //fractional number of neighbors
        double n_r = 0.0;
        for (int ii = 0; ii < neighbors.length; ii++) {
            if (ii == i.label || ii == j.label) {
                this.neighbors[ii] = 1.0;
            } else {
                this.neighbors[ii] = (i.neighbors[ii] + j.neighbors[ii]) / 2;
            }
            n_n += this.neighbors[ii];
            this.remotes[ii] = (i.remotes[ii] + j.remotes[ii]) / 2;
            n_r += this.remotes[ii];
        }
        this.entropy = -(n_n / neighbors.length * lg2(n_n / neighbors.length) + n_r / neighbors.length * lg2(n_r / neighbors.length));

    }

    private double lg2(double d) {
        return Math.log(d) / Math.log(2);
    }
}
